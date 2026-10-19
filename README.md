\# C. elegans Neural Connectome Simulation



An exploratory computational neuroscience project that uses the \*Caenorhabditis elegans\* connectome to simulate neural activity, investigate sensory stimulation, and analyze network response patterns.



\*\*Repository:\*\* https://github.com/sahilc497/Elegans



\## Project Overview



This project loads connectome data through the CECT Python package and implements a rate-based neural network simulation. It investigates how simulated activity propagates when the left or right ASH sensory neuron is stimulated.



The project includes connectome inspection, neural activity simulation, controlled stimulus experiments, left-right response comparisons, parameter robustness tests, coupling ablation, threshold sensitivity analysis, and research reporting.



\## Current Features



\- Connectome data inspection using CECT.

\- Chemical and gap-junction connectivity matrices.

\- Rate-based neural activity simulation.

\- ASHL and ASHR stimulation experiments.

\- Left-right response comparison.

\- Parameter sensitivity and robustness analysis.

\- Coupling ablation experiments.

\- Activity-threshold sensitivity analysis.

\- CSV exports and plots.

\- Consolidated Markdown research report.



\## Technology Stack



\- Python

\- CECT

\- NumPy

\- pandas

\- Matplotlib



\## Project Structure



```text

Elegans/

├── 01\_inspect\_connectome.py

├── 02\_neural\_activity.py

├── ...

├── 09\_connectome\_simulation.py

├── 12\_controlled\_stimuli.py

├── 15\_robustness\_test.py

├── 18\_saturation\_test.py

├── 19\_coupling\_ablation.py

├── 20\_propagation\_analysis.py

├── 21\_threshold\_sensitivity.py

├── 22\_generate\_research\_report.py

├── research\_report/

│   ├── research\_report.md

│   ├── figure\_robustness.png

│   ├── figure\_coupling\_ablation.png

│   └── figure\_threshold\_sensitivity.png

├── requirements.txt

└── .gitignore

```



\## Setup



Python 3.10 is the environment used during development.



\### 1. Clone the repository



```bash

git clone https://github.com/sahilc497/Elegans.git

cd Elegans

```



\### 2. Create a virtual environment



Windows PowerShell:



```powershell

python -m venv .venv

.\\.venv\\Scripts\\Activate.ps1

```



\### 3. Install dependencies



```powershell

python -m pip install --upgrade pip

pip install -r requirements.txt

```



\### 4. Run the simulation



For example:



```powershell

python 09\_connectome\_simulation.py

```



Other numbered scripts perform additional analyses. Run them from the project directory and check their source code for their inputs and output filenames.



\## Research Outputs



The `research\_report/` directory contains the current consolidated analysis, including:



\- Robustness summaries.

\- Coupling ablation results.

\- Threshold sensitivity results.

\- Neural propagation summaries.

\- Figures generated from the simulation outputs.



\## Current Model Limitations



This project is exploratory and should not be interpreted as a validated physiological model of \*C. elegans\*.



\- Chemical connections are treated as positive in the current rate-based model; biological synapse polarity is not comprehensively represented.

\- Gap-junction symmetrization and connectivity normalization are modeling assumptions.

\- Model parameters have not been calibrated against an appropriate experimental neural-activity dataset.

\- Numerical activity thresholds are analytical cutoffs, not measured firing thresholds.

\- Simulated activity patterns do not by themselves establish biological causality or predict worm behavior.



\## Future Work



\- Investigate supported chemical synapse polarity and represent uncertainty explicitly.

\- Validate model outputs against suitable experimental data.

\- Evaluate alternative neuron dynamics and activation functions.

\- Improve experiment reproducibility and automated testing.

\- Extend visualization and comparative analysis.



\## License



No license has been declared yet. Reuse and redistribution permissions should be specified before the project is presented as an openly licensed research resource.



