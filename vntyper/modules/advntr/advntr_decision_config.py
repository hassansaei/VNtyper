"""Pure projection of resolved adVNTR decision and runtime components."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass


@dataclass(frozen=True)
class AdvntrSettings:
    """Typed settings consumed by adVNTR command and result processing."""

    max_frameshift: int
    frameshift_multiplier: int
    vid: int
    output_format: str
    threads: object
    additional_commands: str

    def decision_mapping(self) -> Mapping[str, object]:
        """Return the decision-only settings mapping used by frame processing."""
        return {
            "frameshift_multiplier": self.frameshift_multiplier,
            "max_frameshift": self.max_frameshift,
            "output_format": self.output_format,
            "vid": self.vid,
        }

    def command_mapping(self) -> Mapping[str, object]:
        """Return the combined settings mapping used by command construction."""
        return {
            **self.decision_mapping(),
            "additional_commands": self.additional_commands,
            "threads": self.threads,
        }


def project_advntr_settings(
    decision_component: Mapping[str, object],
    runtime_component: Mapping[str, object],
) -> AdvntrSettings:
    """Project complete resolved components into exact adVNTR stage arguments.

    Args:
        decision_component: Complete ``components.advntr`` mapping.
        runtime_component: Excluded adVNTR runtime mapping.

    Returns:
        Typed combined settings.

    Raises:
        ValueError: If a direct compatibility caller supplies malformed mappings.
    """
    decision = decision_component.get("settings")
    runtime = runtime_component.get("settings")
    if not isinstance(decision, Mapping) or not isinstance(runtime, Mapping):
        raise ValueError("adVNTR decision and runtime settings must be mappings")
    additional_commands = runtime.get("additional_commands")
    if not isinstance(additional_commands, str):
        raise ValueError("adVNTR runtime additional_commands must be a string")
    return AdvntrSettings(
        max_frameshift=int(decision["max_frameshift"]),
        frameshift_multiplier=int(decision["frameshift_multiplier"]),
        vid=int(decision["vid"]),
        output_format=str(decision["output_format"]),
        threads=runtime["threads"],
        additional_commands=additional_commands,
    )
