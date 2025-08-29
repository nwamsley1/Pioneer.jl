#!/bin/bash

# Branch Cleanup Script for Pioneer.jl
# This script helps identify and optionally clean up old branches

echo "=== Git Flow Branch Cleanup Tool ==="
echo ""

# Function to check if branch is merged
is_merged() {
    local branch=$1
    local target=$2
    git merge-base --is-ancestor "$branch" "$target" 2>/dev/null
}

# Fetch latest from remote
echo "Fetching latest from remote..."
git fetch --all --prune

echo ""
echo "=== Local Branches Status ==="
echo ""

# List all local branches and their status
for branch in $(git branch | grep -v "^\*" | sed 's/^ *//'); do
    if [[ "$branch" == "main" || "$branch" == "develop" ]]; then
        continue
    fi
    
    # Check if merged to develop
    if is_merged "$branch" "develop"; then
        echo "✓ $branch (merged to develop)"
    # Check if merged to main
    elif is_merged "$branch" "main"; then
        echo "✓ $branch (merged to main)"
    else
        # Get last commit date
        last_commit=$(git log -1 --format="%cr" "$branch")
        echo "⚠ $branch (unmerged, last commit: $last_commit)"
    fi
done

echo ""
echo "=== Remote Branches Analysis ==="
echo ""

# Analyze remote branches
remote_branches=$(git branch -r | grep -v "HEAD" | sed 's/origin\///' | sort -u)

# Count by category
feature_count=$(echo "$remote_branches" | grep -c "^feature/" || true)
hotfix_count=$(echo "$remote_branches" | grep -c "^hotfix/" || true)
other_count=$(echo "$remote_branches" | grep -v "^feature/\|^hotfix/\|^main$\|^develop$" | wc -l | tr -d ' ')

echo "Feature branches: $feature_count"
echo "Hotfix branches: $hotfix_count"
echo "Other branches: $other_count"

echo ""
echo "=== Branches Not Following Git Flow Convention ==="
echo ""

# List branches not following naming convention
non_gitflow=$(git branch -r | grep -v "HEAD\|main\|develop\|feature/\|hotfix/" | sed 's/.*origin\///')
if [ -n "$non_gitflow" ]; then
    echo "$non_gitflow"
else
    echo "None - all branches follow convention!"
fi

echo ""
echo "=== Recommendations ==="
echo ""

# Check for old branches
old_branches=$(git for-each-ref --format='%(refname:short) %(committerdate:relative)' refs/remotes/origin | \
    grep -v "HEAD\|main\|develop" | \
    awk '$2 == "months" && $1 > 3 {print $1}' | \
    sed 's/origin\///')

if [ -n "$old_branches" ]; then
    echo "Consider deleting these old branches (>3 months):"
    echo "$old_branches"
else
    echo "No branches older than 3 months found."
fi

echo ""
echo "=== Cleanup Commands ==="
echo ""
echo "To delete a merged local branch:"
echo "  git branch -d branch-name"
echo ""
echo "To delete a remote branch:"
echo "  git push origin --delete branch-name"
echo ""
echo "To rename a branch to follow Git Flow:"
echo "  git branch -m old-name feature/new-name"
echo "  git push origin :old-name feature/new-name"
echo ""
echo "To clean up all merged local branches:"
echo "  git branch --merged develop | grep -v 'main\\|develop\\|\\*' | xargs -n 1 git branch -d"