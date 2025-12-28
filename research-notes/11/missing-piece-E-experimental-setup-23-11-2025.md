## 2. Experiment specs for Claude Code (implementation-oriented)

Below is a spec you can give directly to Claude Code. It assumes the usual pattern you already use: an `Experiment` implementation that configures the `Brane`, plus your measurement / diagnostics toolbox.

You can tweak class names and packages, but the logic should be clear.

---

### 2.1 New experiment: `E_Born_ThresholdStatistics`

**Goal**

Implement an experiment that:

1. Initializes a single photon/soliton-like excitation whose envelope splits into two spatial lobes (channels A and B) with configurable relative intensities.
2. Couples two pointer / detector degrees of freedom to disjoint brane regions $\mathcal{R}_A$ and $\mathcal{R}_B$.
3. Uses nonlinear threshold dynamics for the pointers to decide which detector “clicks” first.
4. Repeats the experiment for many microscopic initial conditions that share the same coarse parameters and counts how often A vs B clicks.
5. Optionally scans over different amplitude ratios $|\alpha|^2 : |\beta|^2$.

---

### 2.2 High-level design

**Main pieces:**

1. `BornThresholdExperiment` (implements your experiment interface)
2. `SplitPhotonPhysicalObject` (or reuse your photon object + parameters to create two lobes)
3. `PointerDetector` class to model the pointer DOF with threshold dynamics
4. `BornThresholdConfig` with parameters for amplitudes, detector regions, thresholds, etc.
5. A multi-run harness, e.g. `BornThresholdBatchRunner`, that:

    * runs the experiment many times with different random seeds
    * logs the clicked detector and click time to CSV

---

### 2.3 Configuration class

Create a configuration class for the experiment:

```java
public class BornThresholdConfig {

    // Brane size (elongated in x-direction)
    public int nx = 256;
    public int ny = 32;
    public int nz = 32;

    // Time stepping
    public int totalSteps = 20000;
    public double dt = 1e-3;

    // Photon / soliton parameters
    // amplitudes of the two lobes, alpha^2 + beta^2 should be 1.0
    public double alpha = Math.sqrt(0.5);
    public double beta  = Math.sqrt(0.5);

    // Optional: wavelength, initial position, width, propagation direction
    public double wavelength = /* something at electron/Compton scale in your units */;
    public double packetWidth = /* few grid cells */;
    public double initialX = nx * 0.25;
    public double initialY = ny * 0.5;
    public double initialZ = nz * 0.5;

    // Detector regions (axis-aligned boxes)
    public int detectorAxMin = (int) (nx * 0.05);
    public int detectorAxMax = (int) (nx * 0.15);
    public int detectorBxMin = (int) (nx * 0.85);
    public int detectorBxMax = (int) (nx * 0.95);

    // Use full y,z range for simplicity
    public int detectorAyMin = 0, detectorAyMax = ny - 1;
    public int detectorAzMin = 0, detectorAzMax = nz - 1;
    public int detectorByMin = 0, detectorByMax = ny - 1;
    public int detectorBzMin = 0, detectorBzMax = nz - 1;

    // Pointer parameters
    public double pointerOmega0 = 1.0;   // natural frequency
    public double pointerNonlinear = 1.0; // quartic term strength
    public double pointerDamping  = 0.05;
    public double pointerThreshold = 1.0; // threshold on |Q|

    // Coupling strength between brane and pointer
    public double couplingStrength = 1.0;

    // Noise / microstate parameters
    public double braneNoiseAmplitude = 1e-3;
    public double detectorNoiseAmplitude = 1e-3;

    // Random seed (will be varied between runs)
    public long randomSeed = 12345L;
}
```

Claude Code can extend this as needed.

---

### 2.4 Pointer detector dynamics

Implement a small class for a single pointer with nonlinear threshold dynamics. It does not live on the brane; it is just extra DOF updated in the experiment step.

```java
public class PointerDetector {

    public enum State {
        METASTABLE, CLICKED_POS, CLICKED_NEG
    }

    private double q;      // pointer position
    private double p;      // pointer momentum (conjugate)
    private final double omega0;
    private final double nonlinear;   // coefficient for q^3 term
    private final double damping;
    private final double threshold;   // |q| at which we consider it "clicked"

    private State state = State.METASTABLE;
    private int clickStep = -1;

    public PointerDetector(double omega0, double nonlinear,
                           double damping, double threshold) {
        this.omega0 = omega0;
        this.nonlinear = nonlinear;
        this.damping = damping;
        this.threshold = threshold;
    }

    public void addRandomPerturbation(Random rnd, double amplitude) {
        q += amplitude * (rnd.nextDouble() - 0.5);
        p += amplitude * (rnd.nextDouble() - 0.5);
    }

    /**
     * Advance pointer by one time step, using a simple explicit integrator.
     *
     * drive = coupling term from the brane (e.g. integral over region).
     */
    public void step(double dt, double drive, int stepIndex) {
        if (state != State.METASTABLE) {
            // Once clicked, we no longer evolve it (or we could relax within the basin)
            return;
        }

        // Equation of motion:
        //   dq/dt = p
        //   dp/dt = -omega0^2 * q - nonlinear * q^3 - damping * p + drive
        double dqdt = p;
        double dpdt = - omega0 * omega0 * q
                      - nonlinear * q * q * q
                      - damping * p
                      + drive;

        // Simple Euler or better: semi-implicit / leapfrog if you like
        q += dt * dqdt;
        p += dt * dpdt;

        // Threshold check
        if (Math.abs(q) >= threshold) {
            state = (q > 0) ? State.CLICKED_POS : State.CLICKED_NEG;
            clickStep = stepIndex;
        }
    }

    public boolean isClicked() {
        return state != State.METASTABLE;
    }

    public int getClickStep() {
        return clickStep;
    }

    public State getState() {
        return state;
    }

    public double getQ() {
        return q;
    }

    public double getP() {
        return p;
    }
}
```

You will instantiate one `PointerDetector` for region A and one for region B inside the experiment.

---

### 2.5 Split photon PhysicalObject

You can either:

* Reuse your existing toroidal/photon `PhysicalObject` and pass parameters that produce a directional packet that naturally splits into two lobes, **or**
* Define a simplified `SplitPhotonPhysicalObject` that superposes two Gaussian packets with controlled amplitudes.

Sketch of the simpler version:

```java
public class SplitPhotonPhysicalObject implements PhysicalObject {

    private final BornThresholdConfig cfg;

    public SplitPhotonPhysicalObject(BornThresholdConfig cfg) {
        this.cfg = cfg;
    }

    @Override
    public void initializeBrane(IBrane brane) {
        int nx = brane.getNx();
        int ny = brane.getNy();
        int nz = brane.getNz();

        // Center region where we place the packet
        double x0 = cfg.initialX;
        double y0 = cfg.initialY;
        double z0 = cfg.initialZ;

        double width = cfg.packetWidth;

        // Two lobes: one with positive momentum toward +x (channel B),
        // one with negative momentum toward -x (channel A).
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {

                    double dx = ix - x0;
                    double dy = iy - y0;
                    double dz = iz - z0;
                    double r2 = dx*dx + dy*dy + dz*dz;

                    // Gaussian envelope
                    double envelope = Math.exp(- r2 / (2.0 * width * width));

                    // Assign displacements in X4. Use alpha and beta to
                    // weight left and right propagating parts via phase gradient.
                    // Here we only initialize amplitude; velocity/phase can be set via momentum.
                    double psiA = cfg.alpha * envelope;
                    double psiB = cfg.beta  * envelope;

                    // Combine into one field, but we can encode different momenta via velocity.
                    double xi4 = psiA + psiB;

                    brane.setDisplacement4(ix, iy, iz, xi4);

                    // Set initial velocities in +/- x to create two packets:
                    // e.g. p_x ~ +k for B, -k for A. Implementation depends on your API.
                    double vx = (psiB - psiA) * /* some k-factor */;
                    brane.setVelocity4(ix, iy, iz, vx);
                }
            }
        }
    }
}
```

Claude Code should align this with your actual `IBrane` API (`setDisplacement4`, `setVelocity4` etc. are placeholders; use the real method names).

If you already have a more realistic photon PhysicalObject (with internal torus etc.), you can skip this class and instead adjust its parameters to produce two lobes with relative intensities proportional to `alpha^2` and `beta^2`.

---

### 2.6 The experiment class

Create a new experiment class, e.g. `BornThresholdExperiment`, in the same package as your other experiments:

```java
public class BornThresholdExperiment implements IExperiment {

    private final BornThresholdConfig cfg;
    private final PointerDetector detectorA;
    private final PointerDetector detectorB;
    private final Random rnd;

    private boolean finished = false;
    private String outcome = "NONE"; // "A", "B", or "NONE"

    public BornThresholdExperiment(BornThresholdConfig cfg) {
        this.cfg = cfg;
        this.detectorA = new PointerDetector(
                cfg.pointerOmega0,
                cfg.pointerNonlinear,
                cfg.pointerDamping,
                cfg.pointerThreshold
        );
        this.detectorB = new PointerDetector(
                cfg.pointerOmega0,
                cfg.pointerNonlinear,
                cfg.pointerDamping,
                cfg.pointerThreshold
        );
        this.rnd = new Random(cfg.randomSeed);
    }

    @Override
    public void initializeBraneState(IBrane brane) {
        // 1) create flat brane
        brane.resetToFlatState();

        // 2) initialize split photon / soliton object
        PhysicalObject photon = new SplitPhotonPhysicalObject(cfg);
        photon.initializeBrane(brane);

        // 3) add microscopic noise to brane (hidden variables)
        addBraneNoise(brane, rnd, cfg.braneNoiseAmplitude);

        // 4) add microscopic noise to detectors
        detectorA.addRandomPerturbation(rnd, cfg.detectorNoiseAmplitude);
        detectorB.addRandomPerturbation(rnd, cfg.detectorNoiseAmplitude);
    }

    private void addBraneNoise(IBrane brane, Random rnd, double amplitude) {
        int nx = brane.getNx();
        int ny = brane.getNy();
        int nz = brane.getNz();

        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    double noiseW = amplitude * (rnd.nextDouble() - 0.5);
                    double noiseV = amplitude * (rnd.nextDouble() - 0.5);
                    brane.addDisplacement4(ix, iy, iz, noiseW);
                    brane.addVelocity4(ix, iy, iz, noiseV);
                }
            }
        }
    }

    @Override
    public void step(IBrane brane, int stepIndex) {
        if (finished) return;

        double dt = cfg.dt;

        // 1) compute drive for each detector by integrating over its region
        double driveA = cfg.couplingStrength * computeRegionDrive(
                brane,
                cfg.detectorAxMin, cfg.detectorAxMax,
                cfg.detectorAyMin, cfg.detectorAyMax,
                cfg.detectorAzMin, cfg.detectorAzMax
        );

        double driveB = cfg.couplingStrength * computeRegionDrive(
                brane,
                cfg.detectorBxMin, cfg.detectorBxMax,
                cfg.detectorByMin, cfg.detectorByMax,
                cfg.detectorBzMin, cfg.detectorBzMax
        );

        // 2) advance detectors
        detectorA.step(dt, driveA, stepIndex);
        detectorB.step(dt, driveB, stepIndex);

        // 3) check which detector (if any) has clicked
        boolean aClicked = detectorA.isClicked();
        boolean bClicked = detectorB.isClicked();

        if (aClicked || bClicked) {
            // If both clicked in same step, consider the one with smaller clickStep or higher |Q|
            if (aClicked && !bClicked) {
                outcome = "A";
            } else if (!aClicked && bClicked) {
                outcome = "B";
            } else {
                // tie-breaking rule
                outcome = (Math.abs(detectorA.getQ()) > Math.abs(detectorB.getQ())) ? "A" : "B";
            }
            finished = true;

            // Optionally: record outcome and click times via your MeasurementDevice toolbox
        }

        // Optional: if we exceed totalSteps with no click, mark as "NONE"
        if (stepIndex >= cfg.totalSteps - 1 && !finished) {
            outcome = "NONE";
            finished = true;
        }
    }

    private double computeRegionDrive(IBrane brane,
                                      int xMin, int xMax,
                                      int yMin, int yMax,
                                      int zMin, int zMax) {
        double sum = 0.0;
        int count = 0;

        for (int ix = xMin; ix <= xMax; ix++) {
            for (int iy = yMin; iy <= yMax; iy++) {
                for (int iz = zMin; iz <= zMax; iz++) {
                    // Use displacement in X4 as the drive source; you can
                    // also use local energy density if you have it handy.
                    double w = brane.getDisplacement4(ix, iy, iz);
                    sum += w;
                    count++;
                }
            }
        }

        // Average or just sum; the scaling will be absorbed by cfg.couplingStrength
        return sum / Math.max(1, count);
    }

    public boolean isFinished() {
        return finished;
    }

    public String getOutcome() {
        return outcome;
    }

    public int getClickStepA() {
        return detectorA.getClickStep();
    }

    public int getClickStepB() {
        return detectorB.getClickStep();
    }
}
```

Claude Code will have to adapt method names to your actual `IExperiment` / `IBrane` interfaces, but the structure should be directly portable.

---

### 2.7 Batch runner for statistics

Finally, build a simple batch runner that:

* Varies the random seed (sampling $\lambda_k \in \Gamma_\Psi$),
* Runs the simulation to completion for each seed,
* Logs which detector clicked.

Sketch:

```java
public class BornThresholdBatchRunner {

    private final BornThresholdConfig baseConfig;
    private final int runs;
    private final Path outputCsv;

    public BornThresholdBatchRunner(BornThresholdConfig baseConfig,
                                    int runs,
                                    Path outputCsv) {
        this.baseConfig = baseConfig;
        this.runs = runs;
        this.outputCsv = outputCsv;
    }

    public void runAll() throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(outputCsv, StandardCharsets.UTF_8)) {
            writer.write("runIndex,alpha2,beta2,seed,outcome,clickStepA,clickStepB\n");

            for (int i = 0; i < runs; i++) {
                BornThresholdConfig cfg = cloneConfig(baseConfig);
                cfg.randomSeed = baseConfig.randomSeed + i;

                BornThresholdExperiment experiment = new BornThresholdExperiment(cfg);

                // Create brane with cfg.nx, cfg.ny, cfg.nz and dt
                IBrane brane = BraneFactory.createBrane(cfg.nx, cfg.ny, cfg.nz, cfg.dt);

                // Initialize default brane state (flat) and experiment state
                brane.init(); // or however you currently do it
                experiment.initializeBraneState(brane);

                // Time loop
                for (int step = 0; step < cfg.totalSteps && !experiment.isFinished(); step++) {
                    brane.step();              // update brane by one dt
                    experiment.step(brane, step); // update detectors
                }

                writer.write(String.format(
                        Locale.US,
                        "%d,%.8f,%.8f,%d,%s,%d,%d\n",
                        i,
                        cfg.alpha * cfg.alpha,
                        cfg.beta * cfg.beta,
                        cfg.randomSeed,
                        experiment.getOutcome(),
                        experiment.getClickStepA(),
                        experiment.getClickStepB()
                ));
            }
        }
    }

    private BornThresholdConfig cloneConfig(BornThresholdConfig src) {
        BornThresholdConfig c = new BornThresholdConfig();
        // shallow copy all fields; Claude Code can auto-generate this
        c.nx = src.nx;
        c.ny = src.ny;
        c.nz = src.nz;
        c.totalSteps = src.totalSteps;
        c.dt = src.dt;
        c.alpha = src.alpha;
        c.beta  = src.beta;
        c.wavelength = src.wavelength;
        c.packetWidth = src.packetWidth;
        c.initialX = src.initialX;
        c.initialY = src.initialY;
        c.initialZ = src.initialZ;
        c.detectorAxMin = src.detectorAxMin;
        c.detectorAxMax = src.detectorAxMax;
        c.detectorBxMin = src.detectorBxMin;
        c.detectorBxMax = src.detectorBxMax;
        c.detectorAyMin = src.detectorAyMin;
        c.detectorAyMax = src.detectorAyMax;
        c.detectorAzMin = src.detectorAzMin;
        c.detectorAzMax = src.detectorAzMax;
        c.detectorByMin = src.detectorByMin;
        c.detectorByMax = src.detectorByMax;
        c.detectorBzMin = src.detectorBzMin;
        c.detectorBzMax = src.detectorBzMax;
        c.pointerOmega0 = src.pointerOmega0;
        c.pointerNonlinear = src.pointerNonlinear;
        c.pointerDamping  = src.pointerDamping;
        c.pointerThreshold = src.pointerThreshold;
        c.couplingStrength = src.couplingStrength;
        c.braneNoiseAmplitude = src.braneNoiseAmplitude;
        c.detectorNoiseAmplitude = src.detectorNoiseAmplitude;
        c.randomSeed = src.randomSeed;
        return c;
    }
}
```

To scan different $|\alpha|^2$ values, you can wrap the `BornThresholdBatchRunner` in another loop that modifies `baseConfig.alpha` and `baseConfig.beta = Math.sqrt(1 - alpha^2)`.

---

If you like, in a next step I can help you:

* Tighten the LaTeX section to exactly match your existing notation and labels, and/or
* Align the above Java skeletons with the actual interfaces of your current `Brane`, `PhysicalObject`, and experiment framework once you paste the relevant signatures.
