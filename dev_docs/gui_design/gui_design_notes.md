---
title: "Pioneer Console — visual design notes"
subtitle: "Observations from the running app, 2026-08-12"
geometry: margin=1in
fontsize: 11pt
---

# Scope

These are observations taken from the running SearchDIA page of the Pioneer
Console, cross-checked against the component source. Each item names the
concrete instance that prompted it, so none of this is generic advice.

All items have now been implemented except §8, which was withdrawn after the
claim behind it turned out to be false. Each landed as its own commit so any
one can be judged, or reverted, on its own.

Three measurements in the first draft were taken off Retina screenshots
without dividing by the scale factor, and were wrong by roughly 1.7x. They are
corrected in place below and flagged where they appear.

\newpage

# 1. The number stepper does not fill its field

**Status: done.**

The up/down arrows beside *q-value threshold* leave a white sliver along the
bottom of the field, and the divider between the two arrows sits above the
field's true vertical centre.

This is arithmetic, not a rendering quirk. In `NumField.tsx` the wrapper is

```
display: flex, alignItems: 'stretch', borderRadius: 8, overflow: 'hidden'
```

so the stepper column stretches to the input's full height. But the two
buttons inside it carry a **fixed `height: 18`**. That is 18 + 18 + 1px
divider = 37px of button inside a row whose height is set by the input
(`padding: '9px 11px'` around a 13px mono line, so roughly 41px). The column
is `flexDirection: 'column'` with no `justifyContent`, so the buttons stack at
the top and the leftover falls out of the bottom. The column has no
`background` of its own, so that leftover renders white against the buttons'
`#F5F7F9`.

The fix is to let the buttons divide the column rather than assert a height —
`flex: 1` instead of `height: 18`. That is also more robust: any future change
to the input's font size or padding currently re-opens the gap, whereas
flex-split tracks it automatically.

# 2. One label style per tier

**Status: done.**

Inside *Confidence & Output* two competing label treatments sit side by side
for controls of the same kind:

| control | label style |
|---|---|
| MS data folder (`PathRow`) | 12px, 600, `#344054` |
| Write CSV report (`ToggleRow`) | 13px, 600, `#344054` |
| q-value threshold (`NumField`) | 12px, **regular**, `#475467` |

Both are "a setting with a caption underneath". Rendering the threshold at a
lighter weight and colour makes it read as less important than the toggles
beside it, which inverts the actual hierarchy — the q-value is the more
consequential setting on that card.

The deeper problem was duplication. `LABEL` was defined *identically* in both
`SearchDiaForm.tsx` and `ConvertRawForm.tsx`, while `NumField` carried its own
drifted copy. With three independent definitions there was nothing to stop
them diverging, and they had.

Resolved by extracting `LABEL`, `LABEL_TIGHT` and `HINT` to `lib/styles.ts`,
deleting both duplicates, and pointing all call sites at it. `ToggleRow` also
came down from 13px to 12px — it was the only setting label a point larger
than the rest.

\newpage

# 3. Size fields to their content

The q-value input is 190px wide to hold `0.01`. (An earlier draft said
~320px; that was measured off a Retina capture without dividing by the scale
factor.)

Full-width is correct for the path fields — *MS data folder*, *Spectral
library*, *Results folder* — because paths are long and their length is
unpredictable. But a numeric field that will never hold more than about six
characters reads as unfinished at the same width.

Numerics at 90–120px would look deliberate. It would also shrink the stepper's
visual footprint, so §1 and this compound.

# 4. Let the primary control dominate its row

The `Browse` buttons carry a solid border and bold text, making them visually
heavier than the path inputs they sit beside. The input is the control the
user acts on; Browse is the fallback for when typing a path is inconvenient.

Demoting Browse to a ghost or tertiary treatment would correct the balance
without moving anything.

# 5. Toggle scale

The switches are 38x22 with an 18px knob — not 56x30, as an earlier draft
claimed off an unscaled Retina capture. They still read heavy against 12px
labels, and the dark accent fill against white does as much of that as the
size does.

Trimmed to 34x20 with a 16px knob. A smaller change than the original note
implied, and the fill is the larger factor.

# 6. Two-column blocks need a shared baseline

On *Confidence & Output* the left column ends at `1% FDR = 0.01` while the
right column runs four rows further, leaving a large empty quadrant at the
lower left.

Two ways out: give the left column more to hold, or drop the two-column split
on this card and let the numeric field take its own full-width row above the
toggles. The second is simpler and would pair well with §3.

\newpage

# 7. The title bar is the largest single win

There is a light macOS chrome strip across the top of the window, and the dark
sidebar begins *below* it. The result reads as a web page inside a window
rather than as a native application.

Tauri can extend the sidebar into the title bar — an overlay title-bar style
with the traffic lights inset over the dark chrome. Of everything in this
document, this would do the most for perceived polish, and it is a
configuration change plus a drag region rather than a rework.

Worth noting it interacts with the theme work: the inset traffic lights need
to stay legible on all six palettes, and the light *sand* palette is the one
to check.

# 8. Icon set consistency

**Withdrawn — the premise was wrong.**

The claim was that ConvertRAW used a text arrow while the other two used drawn
icons. It does not: all three are SVGs at matching stroke weight (1.6–1.7).
The `→` in the screenshot is that SVG, and the one in `.raw → .arrow` is
subtitle text.

What survives is much weaker. The arrow spans 14x10 units of its 24 viewBox
while the cylinder spans roughly 14x19 and the magnifier 15x15, so the arrow
occupies less of its box vertically. That is inherent to drawing an arrow
beside two tall closed shapes, and forcing a match would mean a worse arrow.
Left alone deliberately.

# 9. The theme swatch has no affordance

The half-filled circle above *Collapse* gives no signal that it is
interactive. It reads as a stray decoration.

It does already carry a tooltip — the earlier draft missed that — but a
tooltip only helps once you suspect there is something to hover. The word
"theme" beside it, shown at all times rather than only when the picker is
open, is what makes it a control.

\newpage

# Suggested order

Ranked by impact per unit of work:

1. **§7 title bar** — largest visual return, contained change.
2. **§1 stepper** and **§3 field width** — these are what read as "not quite
   right" at a glance.
3. **§2 label hierarchy** — cheap, and it makes the card scan correctly.
4. **§4 Browse weight** and **§5 toggle scale** — refinements once the
   hierarchy above is settled.
5. **§6 column baseline** — partly resolved for free if §3 lands.
6. **§8 icons** and **§9 swatch affordance** — small, independent, do anytime.
