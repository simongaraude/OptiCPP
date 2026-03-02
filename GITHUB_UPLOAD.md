# GitHub Upload Instructions

This document provides step-by-step instructions for uploading the CPP Binder Evolution project to GitHub.

## Prerequisites

1. Git installed on your system
2. GitHub account created
3. Git configured with your credentials

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

## Step 1: Create GitHub Repository

1. Go to https://github.com
2. Click "+" in top right corner
3. Select "New repository"
4. Repository settings:
   - Name: `cpp-binder-evolution`
   - Description: "Evolutionary optimization of cell-penetrating peptide binders"
   - Visibility: Public (or Private)
   - Do NOT initialize with README, .gitignore, or license (we already have these)
5. Click "Create repository"

## Step 2: Prepare Local Repository

```bash
# Navigate to the github directory
cd /path/to/github

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: CPP Binder Evolution v1.0.0

- Production ML CPP predictor (89% accuracy)
- Complete evolutionary algorithm
- Binding hotspot protection
- Structural validation (Step 3.5)
- CPPsite2.0 integration (58 PDB structures)
- Comprehensive documentation"
```

## Step 3: Connect to GitHub

```bash
# Add remote (replace USERNAME with your GitHub username)
git remote add origin https://github.com/USERNAME/cpp-binder-evolution.git

# Verify remote
git remote -v
```

## Step 4: Push to GitHub

```bash
# Push to main branch
git branch -M main
git push -u origin main
```

## Step 5: Verify Upload

1. Go to https://github.com/USERNAME/cpp-binder-evolution
2. Verify all files are present:
   - README.md displaying correctly
   - All directories (binder_evolution/, examples/, pdb_structures/)
   - All documentation files
   - .gitignore working (no __pycache__, .pyc files)

## Step 6: Configure Repository Settings

### Add Topics

Go to repository settings and add topics:
- peptide
- cell-penetrating-peptide
- evolutionary-algorithm
- machine-learning
- drug-design
- protein-engineering
- structural-biology
- bioinformatics

### Add Description

"Evolutionary optimization of cell-penetrating peptide binders with ML guidance and structural validation"

### Add Website

(Optional) If you have documentation hosted elsewhere

### Configure Branch Protection

Settings → Branches → Add rule:
- Branch name pattern: `main`
- Require pull request reviews before merging
- Require status checks to pass before merging

## Step 7: Create Release

### Tag Version

```bash
git tag -a v1.0.0 -m "Release v1.0.0

Production-ready release with:
- ML CPP predictor (89% accuracy)
- Complete evolutionary algorithm
- Structural validation
- CPPsite2.0 integration
- Full documentation"

git push origin v1.0.0
```

### Create GitHub Release

1. Go to repository → Releases → Draft a new release
2. Tag version: v1.0.0
3. Release title: "CPP Binder Evolution v1.0.0"
4. Description:

```markdown
# CPP Binder Evolution v1.0.0

First production-ready release of the CPP Binder Evolution framework.

## Features

- **ML CPP Predictor**: 89% accuracy on 1,075 experimental CPPs
- **Complete Evolution**: All operators (mutation, crossover, selection, elitism)
- **Hotspot Protection**: Preserve binding while optimizing CPP
- **Structural Validation**: Optional Step 3.5 for pose verification
- **CPPsite2.0 Integration**: 58 experimental PDB structures

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```python
from binder_evolution import EvolutionController, create_cpp_scorer

scorer = create_cpp_scorer("cpp_predictor_production.pkl")
controller = EvolutionController(config, scorer, predictor, hotspots)
best, _, _ = controller.run("SQETFSDLWKLLPEN")
```

## Documentation

See README.md for complete documentation.

## Citation

See CITATION.cff for citation information.
```

5. Attach files (optional):
   - cpp_predictor_production.pkl (if not in repository due to size)
   - ml_cpp_predictor.py
   - Example results

6. Click "Publish release"

## Step 8: Add README Badges (Optional)

Edit README.md to add badges at the top:

```markdown
# CPP-Optimized Peptide Binder Evolution

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
```

## Step 9: Set Up CI/CD (Optional but Recommended)

### GitHub Actions for Testing

Create `.github/workflows/tests.yml`:

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, '3.10', 3.11]

    steps:
    - uses: actions/checkout@v2
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v2
      with:
        python-version: ${{ matrix.python-version }}
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install pytest pytest-cov black flake8
    - name: Lint with flake8
      run: |
        flake8 binder_evolution --count --select=E9,F63,F7,F82 --show-source --statistics
    - name: Check code formatting
      run: |
        black --check binder_evolution
    - name: Test with pytest
      run: |
        pytest tests/ --cov=binder_evolution --cov-report=xml
```

## Step 10: Update Documentation

### Create GitHub Wiki (Optional)

1. Go to repository → Wiki → Create first page
2. Add sections:
   - Installation Guide
   - Tutorial
   - API Reference
   - Examples
   - FAQ
   - Troubleshooting

### Create GitHub Pages (Optional)

For hosted documentation:

1. Settings → Pages
2. Source: Deploy from branch
3. Branch: main / docs folder
4. Create docs/ directory with Sphinx or MkDocs

## Step 11: Community Files

GitHub will recognize these files:

- `README.md` - Project overview (✓ included)
- `LICENSE` - License information (✓ included)
- `CONTRIBUTING.md` - Contribution guidelines (✓ included)
- `CODE_OF_CONDUCT.md` - (Create if needed)
- `SECURITY.md` - Security policy (Create if needed)

### Create SECURITY.md (Optional)

```markdown
# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |

## Reporting a Vulnerability

Please report security vulnerabilities to security@example.com

Do not open public issues for security vulnerabilities.
```

## Step 12: Enable Features

### Enable Discussions

Settings → General → Features → Discussions (check)

### Enable Issues

Already enabled by default. Consider creating issue templates:

`.github/ISSUE_TEMPLATE/bug_report.md`:
```markdown
---
name: Bug report
about: Create a report to help us improve
title: '[BUG] '
labels: bug
assignees: ''
---

**Describe the bug**
A clear description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior

**Expected behavior**
What you expected to happen

**Environment:**
- OS: [e.g. Ubuntu 20.04]
- Python version: [e.g. 3.9]
- Package version: [e.g. 1.0.0]
```

## Step 13: Announce Release

1. Update project website (if any)
2. Post on relevant forums:
   - r/bioinformatics
   - r/MachineLearning
   - Biostars
   - ResearchGate
3. Tweet announcement
4. Email collaborators

## Maintenance

### Regular Updates

```bash
# Pull latest changes
git pull origin main

# Create new branch for feature
git checkout -b feature/new-feature

# Make changes, commit
git add .
git commit -m "Add new feature"

# Push branch
git push origin feature/new-feature

# Create pull request on GitHub
```

### Version Updates

When releasing new version:

```bash
# Update version in setup.py
# Update CHANGELOG.md
# Commit changes
git add setup.py CHANGELOG.md
git commit -m "Bump version to 1.1.0"

# Tag new version
git tag -a v1.1.0 -m "Release v1.1.0"
git push origin main --tags

# Create GitHub release
```

## Troubleshooting

### Large Files

If PDB files or model files are too large (>100MB):

```bash
# Install Git LFS
git lfs install

# Track large files
git lfs track "*.pkl"
git lfs track "*.pdb"

# Add .gitattributes
git add .gitattributes
git commit -m "Add Git LFS tracking"
```

### Authentication Issues

For HTTPS:
```bash
# Use personal access token instead of password
git remote set-url origin https://TOKEN@github.com/USERNAME/cpp-binder-evolution.git
```

For SSH:
```bash
# Set up SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"
# Add key to GitHub: Settings → SSH Keys

# Use SSH URL
git remote set-url origin git@github.com:USERNAME/cpp-binder-evolution.git
```

## Final Checklist

- [ ] Repository created on GitHub
- [ ] All files committed and pushed
- [ ] README displays correctly
- [ ] Topics added
- [ ] Release v1.0.0 created
- [ ] License clearly stated
- [ ] Contributing guidelines available
- [ ] Citation file included
- [ ] Documentation complete
- [ ] Examples work
- [ ] .gitignore configured
- [ ] GitHub Actions set up (optional)
- [ ] Issues/Discussions enabled

## Success!

Your repository is now live at:
https://github.com/USERNAME/cpp-binder-evolution

Users can install with:
```bash
git clone https://github.com/USERNAME/cpp-binder-evolution.git
cd cpp-binder-evolution
pip install -r requirements.txt
```

or (after PyPI upload):
```bash
pip install cpp-binder-evolution
```
