#!/bin/bash

# Helper script to automate rebase conflict resolution
# Removes .DS_Store files and continues rebase

while true; do
    echo "=== Rebase Status ==="
    
    # Remove any .DS_Store files that are in conflict
    git rm src/.DS_Store 2>/dev/null
    git rm data/.DS_Store 2>/dev/null
    git rm .DS_Store 2>/dev/null
    
    # Add all changes
    git add -A
    
    # Try to continue rebase
    if git rebase --continue; then
        echo "✅ Rebase completed successfully!"
        break
    fi
    
    # Check if there are conflicts
    if git status --porcelain | grep -q "^UU\|^AA"; then
        echo "⚠️  Non-DS_Store conflicts detected. Manual intervention required."
        echo "Conflicted files:"
        git status --porcelain | grep "^UU\|^AA"
        break
    fi
    
    # Check if rebase is still in progress
    if [ ! -d ".git/rebase-merge" ]; then
        echo "✅ No rebase in progress."
        break
    fi
    
    echo "Continuing to next commit..."
done

echo ""
echo "=== Final Status ==="
git status