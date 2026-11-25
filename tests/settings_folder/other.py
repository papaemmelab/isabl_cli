"""A fixture app that a second module of the folder defines."""

from .base import BaseTool


class OtherTool(BaseTool):
    """Inherit the declarations of `BaseTool`."""
