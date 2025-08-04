#!/bin/bash

# Automated rebase script for Pioneer.jl
# Handles .DS_Store conflicts and preserves feature branch ParameterTuningSearch changes

echo "Starting automated rebase process..."
echo "This will handle .DS_Store conflicts automatically"
echo ""

conflict_count=0
resolved_count=0

while true; do
    # Check if rebase is in progress
    if [ ! -d ".git/rebase-merge" ] && [ ! -d ".git/rebase-apply" ]; then
        echo "✅ No rebase in progress. Done!"
        break
    fi
    
    # Get current rebase status
    if [ -d ".git/rebase-merge" ]; then
        current=$(cat .git/rebase-merge/msgnum 2>/dev/null || echo "?")
        total=$(cat .git/rebase-merge/end 2>/dev/null || echo "?")
        echo "Progress: $current/$total"
    fi
    
    # Check for conflicts
    conflicts=$(git status --porcelain | grep -E "^(UU|AA|DD|DU|UD)")
    
    if [ -n "$conflicts" ]; then
        conflict_count=$((conflict_count + 1))
        echo "Found conflicts (#$conflict_count):"
        echo "$conflicts"
        
        # Handle .DS_Store files
        if echo "$conflicts" | grep -q "\.DS_Store"; then
            echo "  → Removing .DS_Store files..."
            git rm -f **/.DS_Store 2>/dev/null
            git rm -f .DS_Store 2>/dev/null
            git rm -f */.DS_Store 2>/dev/null
            resolved_count=$((resolved_count + 1))
        fi
        
        # Check for ParameterTuningSearch conflicts
        if echo "$conflicts" | grep -q "ParameterTuningSearch"; then
            echo "  → Preserving feature branch version of ParameterTuningSearch files..."
            # Use --theirs to keep feature branch version
            git checkout --theirs src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/*.jl 2>/dev/null
            git add src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/*.jl 2>/dev/null
        fi
        
        # Check for remaining conflicts
        remaining=$(git status --porcelain | grep -E "^(UU|AA)")
        if [ -n "$remaining" ]; then
            echo ""
            echo "⚠️  Non-trivial conflicts remain. Manual intervention required:"
            echo "$remaining"
            echo ""
            echo "Please resolve these conflicts manually, then run this script again."
            break
        fi
    fi
    
    # Add all changes
    git add -A
    
    # Try to continue rebase
    echo "Continuing rebase..."
    if ! git rebase --continue 2>&1 | tee /tmp/rebase_output.txt; then
        # Check if it's just another DS_Store conflict
        if grep -q "\.DS_Store" /tmp/rebase_output.txt; then
            echo "Another .DS_Store conflict detected, handling..."
            continue
        fi
        # Real conflict that needs attention
        echo "Rebase needs attention. Check output above."
    fi
    
    echo ""
done

echo ""
echo "=== Summary ==="
echo "Conflicts encountered: $conflict_count"
echo "Auto-resolved: $resolved_count"
echo ""
echo "Current status:"
git status --short