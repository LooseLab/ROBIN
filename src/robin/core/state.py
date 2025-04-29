"""
Module for storing shared state variables across the application.
"""

from enum import Enum


class ProcessState(Enum):
    """Enum representing the possible states of a process."""

    STARTING = "starting"
    WAITING_FOR_DATA = "waiting_for_data"
    RUNNING = "running"
    STOPPING = "stopping"
    STOPPED = "stopped"


class ProcessType(Enum):
    """Enum representing the possible types of processes."""

    BACKGROUND = "Background"
    PER_FILE = "Per File"
    BATCH = "Batch"
    SYSTEM = "System"


class State:
    """Singleton class to manage shared state across the application."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(State, cls).__new__(cls)
            cls._instance._shutdown_event = False
            cls._instance._process_states = {}  # Dictionary to store process states
            cls._instance._process_types = {}  # Dictionary to store process types
            cls._instance._process_progress = {}  # Dictionary to store process progress
            cls._instance._finished_processes = set()
        return cls._instance

    @property
    def shutdown_event(self):
        return self._shutdown_event

    @shutdown_event.setter
    def shutdown_event(self, value):
        self._shutdown_event = value

    @property
    def process_states(self):
        """Get the dictionary of process states."""
        return self._process_states

    @property
    def process_types(self):
        """Get the dictionary of process types."""
        return self._process_types

    @property
    def process_progress(self):
        """Get the dictionary of process progress values."""
        return self._process_progress

    @property
    def finished_processes(self):
        """Get the set of finished process names."""
        return self._finished_processes

    def start_process(self, process_name: str, process_type: ProcessType):
        """Start tracking a process by name and type."""
        if process_name in self._finished_processes:
            self._finished_processes.remove(process_name)
        self._process_states[process_name] = ProcessState.STARTING
        self._process_types[process_name] = process_type
        self._process_progress[process_name] = 1.0  # Initialize to 100% complete

    def stop_process(self, process_name: str):
        """Stop tracking a process by name and mark it as finished."""
        if process_name in self._process_states:
            del self._process_states[process_name]
            del self._process_types[process_name]
            del self._process_progress[process_name]
            self._finished_processes.add(process_name)

    def set_process_state(self, process_name: str, state: ProcessState):
        """Set the state of a specific process."""
        if process_name in self._process_states:
            self._process_states[process_name] = state

    def get_process_state(self, process_name: str) -> ProcessState:
        """Get the current state of a specific process."""
        return self._process_states.get(process_name, ProcessState.STOPPED)

    def get_process_type(self, process_name: str) -> ProcessType:
        """Get the type of a specific process."""
        return self._process_types.get(process_name)

    def set_process_progress(self, process_name: str, progress: float):
        """Set the progress of a specific process (0.0 to 1.0)."""
        if process_name in self._process_states:
            if not 0.0 <= progress <= 1.0:
                raise ValueError("Progress must be between 0.0 and 1.0")
            self._process_progress[process_name] = progress

    def get_process_progress(self, process_name: str) -> float:
        """Get the progress of a specific process (0.0 to 1.0)."""
        return self._process_progress.get(process_name, 1.0)  # Return 1.0 if not found

    def is_process_running(self, process_name: str) -> bool:
        """Check if a specific process is currently running."""
        return (
            process_name in self._process_states
            and self._process_states[process_name] == ProcessState.RUNNING
        )

    def get_running_process_count(self) -> int:
        """Get the current number of running processes."""
        return sum(
            1
            for state in self._process_states.values()
            if state == ProcessState.RUNNING
        )


# Create a singleton instance
state = State()
