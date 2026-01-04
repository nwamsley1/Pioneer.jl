# Contributing to Pioneer.jl

Thank you for your interest in contributing to Pioneer.jl! This document outlines our development workflow and contribution guidelines.

## Git Flow Model

We follow an adapted version of the [Git Flow](https://nvie.com/posts/a-successful-git-branching-model/) branching model, simplified without release branches.

### Branch Structure

#### Main Branches
- **`main`**: Production-ready code. Only stable, tested code is merged here.
- **`develop`**: Integration branch for features. All development happens here first.

#### Supporting Branches
- **Feature branches** (`feature/*`): For new features and enhancements
- **Hotfix branches** (`hotfix/*`): For critical bug fixes in production

### Development Workflow

#### 1. Working on New Features

```bash
# Start from develop branch
git checkout develop
git pull origin develop

# Create your feature branch
git checkout -b feature/your-feature-name

# Make your changes
# ... code ...

# Commit your changes
git add .
git commit -m "feat: Add your feature description"

# Push to GitHub
git push origin feature/your-feature-name
```

Then create a Pull Request from your feature branch to `develop` on GitHub.

#### 2. Fixing Production Bugs (Hotfixes)

```bash
# Start from main branch for critical fixes
git checkout main
git pull origin main

# Create hotfix branch
git checkout -b hotfix/bug-description

# Fix the bug
# ... code ...

# Commit
git commit -m "fix: Description of the bug fix"

# Push to GitHub
git push origin hotfix/bug-description
```

Create Pull Requests to both `main` and `develop` branches.

### Commit Message Convention

We follow conventional commits format:

- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation changes
- `style:` Code style changes (formatting, etc.)
- `refactor:` Code refactoring
- `test:` Adding or updating tests
- `chore:` Maintenance tasks

Examples:
```
feat: Add support for AlphaPeptDeep model
fix: Correct RT calibration in narrow-window DIA
docs: Update installation instructions
```

### Pull Request Process

1. **Create PR**: Always create PRs to `develop` (or `main` for hotfixes)
2. **Description**: Provide a clear description of changes
3. **Tests**: Ensure all tests pass
4. **Review**: Wait for code review and address feedback
5. **Merge**: PRs will be merged using `--no-ff` to preserve history

### Testing

Before submitting a PR:

```bash
# Run the test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# For specific tests
julia --project=. test/runtests.jl
```

### Documentation

- Update documentation for new features
- Add docstrings to new functions
- Update CLAUDE.md if adding significant new functionality

### Code Style

- Follow Julia style conventions
- Use meaningful variable and function names
- Add comments for complex logic
- Keep functions focused and modular

### Getting Help

- Open an issue for bugs or feature requests
- Use discussions for questions about the codebase
- Check existing issues before creating new ones

### Release Process

Releases are created by maintainers:
1. Merge `develop` → `main` when ready for release
2. Tag the merge commit with version (e.g., `v0.3.0`)
3. GitHub Actions automatically builds and publishes releases


### Branch Naming Examples

Good branch names:
- `feature/add-timstof-support`
- `feature/improve-fdr-calculation`
- `hotfix/fix-memory-leak`
- `hotfix/correct-mass-calculation`

Avoid:
- Generic names like `feature/update` or `fix`
- Personal identifiers like `feature/john-branch`
- Spaces or special characters

## Questions?

Feel free to open a discussion on GitHub or contact the maintainers:
- Nathan Wamsley (wamsleynathan@gmail.com)
- Dennis Goldfarb (dennis.goldfarb@wustl.edu)
