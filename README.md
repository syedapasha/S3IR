<p align="center">
  <h1>S3IR</h1>
  <p align="center">Simulating complex contagion using stochastic simplicial SIR model.</p>
  <p align="center">
    <a href="https://github.com/your-username/S3IR/actions/workflows/ci.yml">
      <img src="https://img.shields.io/badge/build-passing-brightgreen" alt="Build Status">
    </a>
    <a href="LICENSE">
      <img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License">
    </a>
    <a href="https://github.com/your-username/S3IR/pulls">
      <img src="https://img.shields.io/badge/PRs-welcome-brightgreen.svg" alt="PRs Welcome">
    </a>
    <a href="https://github.com/your-username/S3IR/stargazers">
      <img src="https://img.shields.io/github/stars/your-username/S3IR?style=social" alt="GitHub Stars">
    </a>
  </p>
</p>

---

> ### Problem
> Traditional epidemiological models, such as the standard SIR (Susceptible-Infectious-Recovered) framework, often simplify real-world interactions into pairwise relationships. This simplification overlooks the critical role of higher-order social structures in accelerating or mitigating disease transmission. The consequence is often an underestimation or mischaracterization of contagion dynamics, leading to less effective public health interventions and research insights.

> ### Solution
> S3IR (Stochastic Simplicial SIR) transcends these limitations by modeling disease spread on **simplicial complexes**. This approach naturally incorporates multi-person interactions, allowing for a more accurate and nuanced representation of how diseases propagate within complex social networks. By leveraging stochasticity, S3IR provides robust simulations that capture the inherent randomness and variability of real-world contagion, delivering superior predictive power and deeper insights into epidemic dynamics.

---

## Features

*   **Higher-Order Interaction Modeling.** Accurately capture disease transmission through groups and shared contexts, not just pairwise connections, for more realistic simulations.
*   **Stochastic Simulation.** Accommodates the inherent randomness of real-world events, providing a range of probable outcomes and robust statistical analysis for epidemic forecasting.
*   **Comparative Analysis.** Directly compare S3IR's insights against traditional models (e.g., SSEM) to highlight the improved accuracy and nuanced understanding offered by simplicial approaches.
*   **Visual Diagnostics.** Generate insightful visualizations of epidemic progression, network characteristics, and degree distributions to understand complex dynamics at a glance.
*   **Modular MATLAB Design**: Easily extend and adapt the model for custom research scenarios, integrate new datasets, and explore various epidemiological parameters within a powerful numerical environment.
*   **Research-Ready Framework.** Provides a foundation for academic research into complex systems, network science, and epidemiological modeling.

---

## Technical Architecture

| Technology | Purpose                                        | Key Benefit                                                 |
| :--------- | :--------------------------------------------- | :---------------------------------------------------------- |
| MATLAB     | Primary development language for simulations   | Robust mathematical environment, ideal for complex modeling |

### Directory Structure

```
.
├── 📄 Epi_dim.png
├── 📄 deg_dist.m
├── 📄 epi_netw.m
├── 📄 filter_data.m
├── 📄 invs
├── 📄 noEpi_dim.png
├── 📄 s3ir.m
├── 📄 s3ir_vs_ssem.m
├── 📄 sfhh
└── 📄 simp_complex.m
```

---

## Operational Setup

### Prerequisites

To run S3IR, you will need:

*   **MATLAB Runtime Environment**: Version R2017a or newer is recommended for optimal compatibility.

### Installation

Follow these steps to get S3IR up and running on your local machine:

1.  **Clone the Repository**:
    ```bash
    git clone https://github.com/your-username/S3IR.git
    cd S3IR
    ```

2.  **Open in MATLAB**:
    Launch MATLAB and navigate to the `S3IR` directory.

3.  **Run Simulations**:
    Execute the main simulation scripts directly from the MATLAB command window or editor:
    *   To run the core S3IR model:
        ```matlab
        s3ir
        ```
    *   To compare S3IR against the Stochastic Simplicial Epidemic Model (SSEM):
        ```matlab
        s3ir_vs_ssem
        ```

### Environment Configuration

This project does not require external `.env` files or specific configuration files beyond standard MATLAB setup. All parameters are managed directly within the `.m` scripts.

---

