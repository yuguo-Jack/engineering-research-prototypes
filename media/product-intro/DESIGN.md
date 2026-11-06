# Engineering Research Prototypes — video identity

## Style Prompt

V2 is an exactly 30-second product launch trailer, cut to 140 BPM. Keep the
website's warm white, deep green, teal and orange, but alternate warm-white
evidence scenes with deep-green kinetic scenes. Nine chapters contain twelve
visual beats: a 1.7-second hook, a snap reveal, diagram builds, a 3D document
stack, domain line art, two real CFD cases and an invitation. Use perspective,
masked whip cuts, impulse lines, synchronized strokes and deliberate push-ins.
Create energy through motion and rhythm instead of filling the frame with text.
Use real project outputs for both CFD examples and label future domains as
directions, rather than suggesting that unbuilt prototypes already exist.

## Colors

- `#fbfcfb` — the website's primary light canvas.
- `#172d30` — primary ink, headline and diagram labels.
- `#0b746b` — main accent, active paths, key figures and selected states.
- `#ca713d` — small signal accents and key milestones.
- `#dce5e2` — structural rules and quiet diagram background.

Secondary text may use the website's `#637477`; any derived tints must remain
in this palette. Evidence remains on light paper; hooks, transitions and the
closing use the existing deep-green ink as their background.

## Typography

- Chinese and interface text: `ERP Sans`, an embedded and renamed subset of
  the OFL-licensed Noto Sans SC variable font. Its copyright and license are
  retained in `fonts/OFL.txt`; no proprietary font files are distributed.
- Numbers, paths and English labels: `Segoe UI` / `Consolas`, matching the site.
- 1920 × 1080; short titles 72–116 px and numbers 104–156 px. Keep reading
  copy concise. Text and illustrations occupy separate columns or bands.
- Main picture stays above y=820; a dedicated caption rail at y=884–960
  contains short 30–36 px captions. Never place subtitles over figures.
- Native player captions are optional and off by default because the film has
  its own caption rail. Detail and reproducibility conditions remain on the site.

## Motion

Use GSAP with deterministic paused timelines. Entrance choreography should mix
three compatible eases per scene and combine translation with opacity/scale.
Transitions take 0.18–0.26 seconds: diagonal masks, hard cuts on musical beats
and compact pushes. Use actual 3D transforms on document planes and rotating
mesh line art. Build key figures in time with drum hits. Avoid 5-second holds.
Transitions handle scene exits; do not pre-fade scenes. Each case has one
dominant graphic and at most three concise facts, with no long overlay title.

## What NOT to Do

- No purple/cyan neon, generic dark-tech gradients or unrelated stock footage.
- No fabricated solver results, official scores or same-machine speedup claims.
- Do not present prospective radar/robotics/embedded directions as completed work.
- Do not replace M6's confirmed window time with its earlier starting timestamp.
- No donation promotion in the product narrative; the closing invites real needs.
- No rendered caches, intermediate frames or temporary audio in Git.
