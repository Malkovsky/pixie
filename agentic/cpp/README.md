# Shared C++ Agent Skills

This subtree contains reusable C++ agent skills and related commands.

Keep this tree generic:

- Do not add project-specific benchmark names, CMake options, or paths.
- Keep reusable scripts beside the skills that use them.
- Put project-specific examples in the consuming repository under
  `agentic/local/cpp/skills/<skill-name>/EXAMPLES.md`.

When using a skill in a project, read:

1. `agentic/cpp/skills/<skill-name>/SKILL.md`
2. `agentic/local/cpp/skills/<skill-name>/EXAMPLES.md`, if present
