# Interactive Examples

## Overview

This section contains Jupyter notebooks and Python scripts demonstrating practical QuREX applications. These provide hands-on experience with real data and complete workflows.

## Available Notebooks

### [Exciton-Phonon Coupling Analysis](exciton_phonon)
**Type**: Jupyter Notebook | **Level**: Intermediate

Analysis of electron-phonon interactions in 2D materials:
- Workflow from DFT to exciton-phonon matrix elements
- Temperature-dependent optical properties
- Phonon-assisted absorption and emission processes
- Real-space visualization of coupling mechanisms

**Materials**: MoS₂, hBN, graphene  
**Prerequisites**: Basic Python, GW+BSE calculations

### [Exciton Group Theory Classification](exciton_group_theory_example)
**Type**: Python Script | **Level**: Advanced

Symmetry analysis and classification of excitonic states:
- Point group identification and character tables
- Irreducible representation decomposition
- Optical selection rules determination
- Bright/dark exciton classification

**Prerequisites**: Group theory basics, BSE calculations with symmetry data

## Getting Started

### Environment Setup

```bash
# Install dependencies
pip install jupyter matplotlib numpy scipy yambopy

# Launch Jupyter
jupyter notebook
```

### Data Requirements

Each notebook includes:
- Sample data for immediate use
- Complete input files
- Reference results for validation
- Download instructions for additional datasets

## Notebook Structure

### 1. Introduction
- Theoretical background and key concepts
- Physical context and learning objectives

### 2. Setup
- Library imports and dependencies
- Data loading and parameter definition

### 3. Analysis
- Step-by-step calculations
- Visualization with explanations
- Physical interpretation of results

### 4. Extensions
- Additional analysis possibilities
- Customization for different systems
- Research applications

## Usage Patterns

### Learning Mode
- Run cells sequentially to understand workflow
- Modify parameters to explore effects
- Experiment with different materials

### Research Mode
- Adapt notebooks for your own data
- Extend analysis techniques
- Generate publication figures

### Reference Mode
- Find specific analysis techniques
- Extract useful code snippets
- Learn best practices

## Technical Support

### Common Issues
- **Installation**: Dependency conflicts and solutions
- **Data Loading**: File format and path issues
- **Calculations**: Debugging computational problems
- **Visualization**: Plot rendering issues

### Troubleshooting

Check environment:
```python
import sys, numpy, matplotlib, scipy
print(f"Python version: {sys.version}")
print("Core packages loaded successfully")
```

Validate data files:
```python
import os
required_files = ['ndb.BS_diago_Q1', 'ns.db1', 'ndb.elph']
for file in required_files:
    print(f"{'✓' if os.path.exists(file) else '✗'} {file}")
```

## Contributing

We welcome community contributions:
- New notebooks for different applications
- Bug fixes and improvements
- Enhanced documentation
- Performance optimizations

Process:
1. Fork repository
2. Develop content
3. Test thoroughly
4. Submit pull request