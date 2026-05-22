---
name: Feature Request (Developer)
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''
type: 'Feature'
---

<!--
**For internal developers only.** External contributors should use the external issue template.
Please do not post usage questions here. Ask them on Discord or contact tyndp@openenergytransition.org.
-->

## Parent Issue (if sub-issue)
<!-- Is this a sub-issue of a larger parent issue? If yes, link the parent here and **edit the title above to start with [SUB]**. -->
<!-- Example: Relates to #123 -->
<!-- Leave blank if this is a standalone feature request -->


## Description
<!-- Describe the requested feature or change. Extend the issue title with context if needed. -->


## What (Expected Impact)
<!-- What problem does this solve, or what value does it add? Be brief. -->


## Expected Impact
<!-- What will change? Systems, APIs, data flows, dependencies affected (if relevant). Mention affected files, workflows, datasets, scenarios, or outputs if known. Also feel free to use the labels available (such as SB, CBA, documentation, etc). -->

## Change Classification
<!-- Pick one. This determines the approval path according to OET-ISMS-006, Annex A. -->

- [ ] **Standard** — Low-risk, routine (e.g. docs update, non-critical dependency bump) → *Developer self-approves*
- [ ] **Normal** — Logic, infra, or API changes (e.g. new endpoint, schema change) → *Project Lead approval before work starts*
- [ ] **Emergency (Break-Glass)** — Critical hotfix / active security breach → *Verbal approval from Senior Lead; ticket updated retrospectively within 24h*


## Security & Data Risk Check
<!-- Please answer honestly. Checking a box does not block the change; it flags what needs review.
Leave unchecked if the answer is No or not applicable. -->

- [ ] This touches sensitive data (user PII, client grid data, credentials).
- [ ] This adds or changes external data sources, download URLs, APIs, services, credentials, or tokens.
- [ ] This changes the attack surface (new ports, APIs, external dep, auth flows).


**Risk level** (required for major features — skip for minor changes):
<!-- High risk requires a kick-off security review. If personal data is involved, assess whether a PIA is required. -->
- [ ] Low
- [ ] Medium
- [ ] High

## Tasks
<!-- List any specific tasks that need to be done for this task (if known). -->

- [ ] 
- [ ] 
- [ ] 

## Test & Deployment Notes


## Test & Validation Plan
<!-- Describe how the change(s) should be tested or checked. Can list concrete commands to run, files to check, expected outputs, etc. -->
<!-- Examples: relevant pixi/snakemake commands, benchmark comparisons, etc. -->

## Rollback / Recovery Plan
<!-- Required for Normal and Emergency changes. Briefly describe how to undo the change within about 15 minutes if it breaks something. -->
<!-- Example: revert PR commit; restore previous config; pin previous dataset version; disable feature flag. -->

## Additional Context
<!-- Links, references, implementation ideas, screenshots, related issues/PRs, etc. -->