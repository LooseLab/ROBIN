"""
Module for storing shared state variables across the application.
"""

class State:
    """Singleton class to manage shared state across the application."""
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(State, cls).__new__(cls)
            cls._instance._shutdown_event = False
            cls._instance._running_processes = set()
            cls._instance._finished_processes = set()
        return cls._instance
    
    @property
    def shutdown_event(self):
        return self._shutdown_event
        
    @shutdown_event.setter
    def shutdown_event(self, value):
        self._shutdown_event = value

    @property
    def running_processes(self):
        """Get the set of currently running process names."""
        return self._running_processes

    @property
    def finished_processes(self):
        """Get the set of finished process names."""
        return self._finished_processes

    def start_process(self, process_name: str):
        """Start tracking a process by name."""
        if process_name in self._finished_processes:
            self._finished_processes.remove(process_name)
        self._running_processes.add(process_name)

    def stop_process(self, process_name: str):
        """Stop tracking a process by name and mark it as finished."""
        if process_name in self._running_processes:
            self._running_processes.remove(process_name)
            self._finished_processes.add(process_name)

    def is_process_running(self, process_name: str) -> bool:
        """Check if a specific process is currently running."""
        return process_name in self._running_processes

    def get_running_process_count(self) -> int:
        """Get the current number of running processes."""
        return len(self._running_processes)

# Create a singleton instance
state = State() 