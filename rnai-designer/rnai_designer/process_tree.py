"""Contain a worker and descendants; no scientific timeout policy here."""
import os
import signal
import time


class ProcessTree:
    def __init__(self):
        self.process = None
        self.job = None
        if os.name == 'nt':
            self.job = WindowsJob()

    def attach(self, process):
        # The worker waits for stdin before starting any work/descendants.
        if self.job:
            self.job.assign(process)
        self.process = process

    def stop(self):
        if self.job:
            self.job.terminate()
        elif self.process:
            try:
                os.killpg(self.process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        if self.process:
            self.process.wait(timeout=5)
        if self.job:
            until = time.monotonic()+5
            while self.job.active():
                if time.monotonic() >= until:
                    raise OSError('No se pudo confirmar la terminacion del arbol de procesos')
                time.sleep(.02)

    def close(self):
        if self.job:
            self.job.close()


class WindowsJob:
    """Windows 8+ Job Object, no breakaway, kill on last handle close."""
    def __init__(self):
        import ctypes as c
        from ctypes import wintypes as w
        class BasicLimit(c.Structure):
            _fields_ = [('PerProcessUserTimeLimit', c.c_longlong), ('PerJobUserTimeLimit', c.c_longlong),
                        ('LimitFlags', w.DWORD), ('MinimumWorkingSetSize', c.c_size_t),
                        ('MaximumWorkingSetSize', c.c_size_t), ('ActiveProcessLimit', w.DWORD),
                        ('Affinity', c.c_size_t), ('PriorityClass', w.DWORD), ('SchedulingClass', w.DWORD)]
        class IO(c.Structure):
            _fields_ = [(key,c.c_ulonglong) for key in ('ReadOperationCount','WriteOperationCount',
                        'OtherOperationCount','ReadTransferCount','WriteTransferCount','OtherTransferCount')]
        class Extended(c.Structure):
            _fields_ = [('BasicLimitInformation', BasicLimit), ('IoInfo', IO),
                        ('ProcessMemoryLimit',c.c_size_t),('JobMemoryLimit',c.c_size_t),
                        ('PeakProcessMemoryUsed',c.c_size_t),('PeakJobMemoryUsed',c.c_size_t)]
        class Accounting(c.Structure):
            _fields_ = [(key,c.c_longlong) for key in ('TotalUserTime','TotalKernelTime',
                          'ThisPeriodTotalUserTime','ThisPeriodTotalKernelTime')] + [
                        (key,w.DWORD) for key in ('TotalPageFaultCount','TotalProcesses',
                          'ActiveProcesses','TotalTerminatedProcesses')]
        self.c, self.Accounting = c, Accounting
        self.api = c.WinDLL('kernel32',use_last_error=True)
        specs = {'CreateJobObjectW': ([c.c_void_p,w.LPCWSTR],w.HANDLE),
                 'SetInformationJobObject': ([w.HANDLE,c.c_int,c.c_void_p,w.DWORD],w.BOOL),
                 'AssignProcessToJobObject': ([w.HANDLE,w.HANDLE],w.BOOL),
                 'TerminateJobObject': ([w.HANDLE,w.UINT],w.BOOL),
                 'QueryInformationJobObject': ([w.HANDLE,c.c_int,c.c_void_p,w.DWORD,c.c_void_p],w.BOOL),
                 'CloseHandle': ([w.HANDLE],w.BOOL)}
        for name,(args,result) in specs.items():
            fn=getattr(self.api,name);fn.argtypes=args;fn.restype=result
        self.handle=self.api.CreateJobObjectW(None,None)
        if not self.handle:
            raise c.WinError(c.get_last_error())
        limits=Extended();limits.BasicLimitInformation.LimitFlags=0x2000  # KILL_ON_JOB_CLOSE
        if not self.api.SetInformationJobObject(self.handle,9,c.byref(limits),c.sizeof(limits)):
            error=c.WinError(c.get_last_error());self.close();raise error

    def assign(self, process):
        if not self.api.AssignProcessToJobObject(self.handle,int(process._handle)):
            raise self.c.WinError(self.c.get_last_error())

    def active(self):
        info=self.Accounting()
        if not self.api.QueryInformationJobObject(self.handle,1,self.c.byref(info),self.c.sizeof(info),None):
            raise self.c.WinError(self.c.get_last_error())
        return info.ActiveProcesses

    def terminate(self):
        if not self.api.TerminateJobObject(self.handle,124):
            raise self.c.WinError(self.c.get_last_error())

    def close(self):
        if self.handle:
            self.api.CloseHandle(self.handle)
            self.handle=None
