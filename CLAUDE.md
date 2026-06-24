# CLAUDE.md

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

## 5. Julia Performance & Memory Discipline

**Pioneer is performance- and allocation-critical (per-scan/per-precursor hot loops, threaded deconvolution). When editing Pioneer code, follow Julia's official guidance:**
- Performance Tips: https://docs.julialang.org/en/v1/manual/performance-tips/
- Memory Management: https://docs.julialang.org/en/v1/manual/memory-management/

Apply these *especially in hot paths* (anything inside per-scan, per-precursor, or per-fragment loops). Outside hot paths, prefer clarity (§2) — don't micro-optimize where the profiler doesn't point.

**Type stability**
- Functions should return a single concrete type; avoid code where a variable's type changes. Verify hot functions with `@code_warntype` (no red `Any`/`Union`).
- Use function barriers to isolate unavoidable instability (e.g. reading untyped data) from the hot kernel.

**Concrete types**
- No abstract field types in structs and no abstract element types in containers (`Vector{AbstractFoo}` dispatches dynamically and boxes). Use concrete or parametric types.
- Module-level globals used in hot code must be `const` (or a `Ref`/typed constant).

**Allocation (the usual win here)**
- Don't allocate in hot loops. Pre-allocate buffers/workspaces once and reuse them (per-thread where threaded); prefer in-place `!` functions.
- Use `@view`/`view(...)` instead of array slices (`a[1:n]` copies).
- Small fixed-size collections (≤ ~16) → `NTuple` or `StaticArrays.MVector` (stack, alloc-free), not heap `Vector`. Mark the producing function `@inline` so non-escaping stack scratch elides.
- Know stack vs heap: `isbits`/immutable values are stack-allocated and allocation-free; mutable/heap types allocate and add GC pressure.

**Measure, don't guess**
- Use `@time`/`@allocated`/`@btime`. Distinguish *cumulative* allocation (`@allocated`, `@timed.bytes` — process-global, sums across threads) from *peak working set* (RSS): a large cumulative number against a small RSS means avoidable churn, not a large footprint.
- `@allocated` on a standalone call is pessimistic (forces a non-inlined boundary); measure the realistic call pattern (e.g. amortized over a loop) before concluding a function still allocates.
- Profile first; optimize where the data points. State the before/after allocation or timing number when claiming a perf win.

**Safety**
- Only add `@inbounds`/`@simd` when bounds are provably safe, and verify results are bit-identical to the unoptimized version.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.