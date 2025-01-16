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

### Reporting Bugs

If you find a bug in **VNtyper**, please help us by creating a new issue:

1. **Check for Existing Issues**: Before submitting a new bug report, please search the [Issues](https://github.com/hassansaei/VNtyper/issues) to see if it has already been reported.
2. **Open a New Issue**: If you don’t find an existing report, click on the `New Issue` button.
3. **Use a Clear Title**: Provide a short, descriptive title of the bug.
4. **Describe the Bug**:
   - Explain the issue clearly and concisely.
   - Include steps to reproduce the bug.
   - Mention any relevant configuration or environment details (e.g., OS, Python version).
   - Describe the expected behavior vs. the actual behavior.
5. **Add Additional Context (if necessary)**: Screenshots, logs, or error messages can be very helpful.

By reporting bugs promptly and clearly, you help us investigate and resolve them more efficiently.

### Suggesting Enhancements

We welcome ideas to improve **VNtyper**! To suggest a new feature or enhancement:

1. **Check for Existing Requests**: Search the [Issues](https://github.com/hassansaei/VNtyper/issues) for similar feature requests.
2. **Open a New Issue**: If the idea is not already discussed, click `New Issue`.
3. **Use a Descriptive Title**: Clearly state what feature or improvement you would like to see.
4. **Provide Details**:
   - Explain why the enhancement would be useful.
   - Suggest any potential approach, relevant references, or examples.
5. **Discuss and Refine**: Engage with the community and maintainers on the issue to refine the idea and agree on an approach.

### Submitting Pull Requests

When you're ready to contribute code:

1. **Fork the Repository**: Click on the `Fork` button at the top right of the repository page.

2. **Clone Your Fork**:

   ```bash
   git clone https://github.com/your-username/vntyper.git
   cd vntyper
   ```

3. **Create a New Branch**:

   ```bash
   git checkout -b feature/your-feature-name
   ```

4. **Make Changes**: Implement your feature or fix.

5. **Commit Changes**: Follow the [Commit Messages](#commit-messages) guidelines.

6. **Push to Your Fork**:

   ```bash
   git push origin feature/your-feature-name
   ```

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

   ```bash
   git clone https://github.com/your-username/vntyper.git
   cd vntyper
   ```

2. **Install Mamba (if not already installed)**:

   ```bash
   conda install mamba -n base -c conda-forge
   ```

3. **Create a Development Environment**:

   ```bash
   mamba env create -f environment.yml
   conda activate vntyper-dev
   ```

   *Note*: Ensure that an `environment.yml` file exists with the necessary development dependencies.

4. **Install VNtyper in Editable Mode**:

   ```bash
   pip install -e .
   ```

---

## Coding Guidelines

### Code Style

We follow the **PEP 8** style guide with some project-specific conventions.

- **Use `black`** for code formatting:

  ```bash
  black vntyper/
  ```

- **Lint Your Code** with `flake8`:

  ```bash
  flake8 vntyper/
  ```

- **Type Checking** with `mypy` (if applicable):

  ```bash
  mypy vntyper/
  ```

### Commit Messages

- **Format**: Use the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) style.
- **Structure**:

  ```
  type(scope): subject

  body (optional)

  footer (optional)
  ```

- **Example**:

  ```
  feat(pipeline): add new genotyping method

  Implemented a new method for genotyping using advanced algorithms.

  Closes #123
  ```

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