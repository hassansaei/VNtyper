# Contributing to VNtyper

Thank you for considering contributing to **VNtyper**! Your contributions help make this project better for everyone. This document outlines the process for contributing, reporting issues, setting up the development environment, and coding standards.

## Table of Contents

1. [How Can I Contribute?](#how-can-i-contribute)
   - [Reporting Bugs](#reporting-bugs)
   - [Suggesting Enhancements](#suggesting-enhancements)
   - [Submitting Pull Requests](#submitting-pull-requests)
2. [Development Setup](#development-setup)
   - [Prerequisites](#prerequisites)
   - [Setting Up the Development Environment](#setting-up-the-development-environment)
3. [Coding Guidelines](#coding-guidelines)
   - [Code Style](#code-style)
   - [Commit Messages](#commit-messages)
4. [Pull Request Process](#pull-request-process)
5. [Community Guidelines](#community-guidelines)

---

## How Can I Contribute?

### Submitting Pull Requests

When you're ready to contribute code:

1. **Fork the Repository**: Click on the `Fork` button at the top right of the repository page.

2. **Clone Your Fork**:

   ~~~bash
   git clone https://github.com/your-username/vntyper.git
   cd vntyper
   ~~~

3. **Create a New Branch**:

   ~~~bash
   git checkout -b feature/your-feature-name
   ~~~

4. **Make Changes**: Implement your feature or fix.

5. **Commit Changes**: Follow the [Commit Messages](#commit-messages) guidelines.

6. **Push to Your Fork**:

   ~~~bash
   git push origin feature/your-feature-name
   ~~~

7. **Open a Pull Request**: Go to the original repository and click on `New Pull Request`.
   - Provide a clear description of your changes.
   - Reference any related issues (`Closes #issue-number`).

---

## Development Setup

### Prerequisites

- **Python**: Version 3.9 or higher
- **Git**: Version control system
- **Conda/Mamba**: For environment management

### Setting Up the Development Environment

1. **Clone the Repository**:

   ~~~bash
   git clone https://github.com/your-username/vntyper.git
   cd vntyper
   ~~~

2. **Install Mamba (if not already installed)**:

   ~~~bash
   conda install mamba -n base -c conda-forge
   ~~~

3. **Create a Development Environment**:

   ~~~bash
   mamba env create -f environment.yml
   conda activate vntyper-dev
   ~~~

   *Note*: Ensure that an `environment.yml` file exists with the necessary development dependencies.

4. **Install VNtyper in Editable Mode**:

   ~~~bash
   pip install -e .
   ~~~

---

## Coding Guidelines

### Code Style

We follow the **PEP 8** style guide with some project-specific conventions.

- **Use `black`** for code formatting:

  ~~~bash
  black vntyper/
  ~~~

- **Lint Your Code** with `flake8`:

  ~~~bash
  flake8 vntyper/
  ~~~

- **Type Checking** with `mypy` (if applicable):

  ~~~bash
  mypy vntyper/
  ~~~

### Commit Messages

- **Format**: Use the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) style.
- **Structure**:

  ~~~
  type(scope): subject

  body (optional)

  footer (optional)
  ~~~

- **Example**:

  ~~~
  feat(pipeline): add new genotyping method

  Implemented a new method for genotyping using advanced algorithms.

  Closes #123
  ~~~

---

## Pull Request Process

1. **Ensure All Tests Pass**: Before submitting, make sure all tests pass locally.
2. **Update Documentation**: If your changes affect documentation, update it accordingly.
3. **Follow the Template**: Use the provided pull request template and fill out all sections.
4. **Address Review Comments**: Be responsive to feedback and make necessary changes.

---

## Community Guidelines

- Be respectful and considerate.
- Provide constructive feedback.
- Follow the project's [Code of Conduct](CODE_OF_CONDUCT.md) (if available).

---
