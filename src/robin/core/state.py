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

class State:
    """Singleton class to manage shared state across the application."""
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(State, cls).__new__(cls)
            cls._instance._shutdown_event = False
            cls._instance._process_states = {}  # Dictionary to store process states
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
    def finished_processes(self):
        """Get the set of finished process names."""
        return self._finished_processes

    def start_process(self, process_name: str):
        """Start tracking a process by name."""
        if process_name in self._finished_processes:
            self._finished_processes.remove(process_name)
        self._process_states[process_name] = ProcessState.STARTING

    def stop_process(self, process_name: str):
        """Stop tracking a process by name and mark it as finished."""
        if process_name in self._process_states:
            del self._process_states[process_name]
            self._finished_processes.add(process_name)

    def set_process_state(self, process_name: str, state: ProcessState):
        """Set the state of a specific process."""
        if process_name in self._process_states:
            self._process_states[process_name] = state

    def get_process_state(self, process_name: str) -> ProcessState:
        """Get the current state of a specific process."""
        return self._process_states.get(process_name, ProcessState.STOPPED)

    def is_process_running(self, process_name: str) -> bool:
        """Check if a specific process is currently running."""
        return process_name in self._process_states and self._process_states[process_name] == ProcessState.RUNNING

    def get_running_process_count(self) -> int:
        """Get the current number of running processes."""
        return sum(1 for state in self._process_states.values() if state == ProcessState.RUNNING)

# Create a singleton instance
state = State() 