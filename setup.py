from setuptools import setup, find_packages

setup(
    name="branesim",
    version="0.1.0",
    description="PyTorch-based Brane simulation with Velocity Verlet integration",
    author="Lukas Molzberger",
    packages=find_packages(),
    install_requires=[
        "torch>=2.0.0",
        "numpy>=1.24.0",
        "matplotlib>=3.7.0",
        "scipy>=1.10.0",
    ],
    python_requires=">=3.8",
)
