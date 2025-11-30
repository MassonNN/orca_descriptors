# GitHub Workflows

This directory contains GitHub Actions workflows for the ORCA Descriptors project.

## Workflows

### 1. Build Documentation (`docs.yml`)

Builds and deploys documentation to GitHub Pages.

**Triggers:**
- Push to `main`/`master` branch (when `docs/` changes)
- Pull requests to `main`/`master` branch
- Manual trigger via `workflow_dispatch`

**What it does:**
- Builds English and Russian documentation using Sphinx
- Combines both language versions into a single structure
- Deploys to GitHub Pages
- Documentation structure:
  - English: `/` (root)
  - Russian: `/ru/`

**Requirements:**
- GitHub Pages must be enabled in repository settings
- Set source to "GitHub Actions" in Pages settings (Settings → Pages → Source: GitHub Actions)

**Access:**
- Documentation will be available at: `https://<username>.github.io/orca_descriptors/`
- English: `https://<username>.github.io/orca_descriptors/`
- Russian: `https://<username>.github.io/orca_descriptors/ru/`

### 2. Publish to PyPI (`publish.yml`)

Publishes the package to PyPI using OpenID Connect (OIDC) authentication (trusted publishing).

**Triggers:**
- Release published (when a new release is created on GitHub)
- Manual trigger via `workflow_dispatch`

**What it does:**
- Builds the package using Poetry
- Publishes to PyPI using trusted publishing (OIDC)
- No API tokens required when using trusted publishing

**Setup Instructions:**

#### Step 1: Enable Trusted Publishing on PyPI

1. Go to https://pypi.org/manage/account/publishing/
2. Click "Add a new pending publisher"
3. Fill in the form:
   - **PyPI project name**: `orca-descriptors` (must match your package name in `pyproject.toml`)
   - **Owner**: Your GitHub username or organization name
   - **Repository name**: `orca_descriptors` (or your actual repository name)
   - **Workflow filename**: `.github/workflows/publish.yml`
   - **Environment name (optional)**: `pypi` (recommended for security)
4. Click "Add"
5. The publisher will be in "pending" state until you run the workflow for the first time

#### Step 2: Create GitHub Environment (Recommended)

1. Go to your repository: Settings → Environments
2. Click "New environment"
3. Name it `pypi` (must match the environment name in the workflow)
4. Optionally add protection rules:
   - Required reviewers (for production safety)
   - Deployment branches (restrict to specific branches)
5. Click "Configure environment"

#### Step 3: Publishing

**Automatic (Recommended):**
1. Update version in `pyproject.toml`
2. Commit and push changes
3. Create a new release on GitHub:
   - Go to Releases → "Draft a new release"
   - Create a new tag (e.g., `v0.1.0`)
   - Publish the release
4. The workflow will automatically trigger and publish to PyPI

**Manual:**
1. Go to Actions → "Publish to PyPI"
2. Click "Run workflow"
3. Select branch (usually `main` or `master`)
4. Enter version number (must match `pyproject.toml`)
5. Click "Run workflow"

**Note:** The workflow uses OIDC (OpenID Connect) for authentication via trusted publishing. This is more secure than API tokens and doesn't require storing secrets. The `trusted-publishing: true` option enables this feature.

## Troubleshooting

### Documentation not deploying

- Check that GitHub Pages is enabled: Settings → Pages → Source should be "GitHub Actions"
- Verify the workflow has `pages: write` permission
- Check workflow logs for errors

### PyPI publishing fails

- Verify trusted publishing is configured correctly on PyPI
- Check that the environment name matches (`pypi`)
- Ensure the PyPI project name matches `pyproject.toml` (`orca-descriptors`)
- Verify the repository name and owner are correct
- Check workflow logs for detailed error messages

### Version mismatch

- Ensure the version in `pyproject.toml` matches the release tag
- For manual publishing, enter the exact version from `pyproject.toml`
