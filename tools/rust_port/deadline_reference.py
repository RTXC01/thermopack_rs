"""An explicit execution deadline for original solvers that may never return.

Timeouts are censored observations, not proof of a nonconvergent source result.
The worker is intentionally killed and reaped; the case is never retried.
"""
from queue import Queue,Empty
from threading import Thread
from reference import ReferenceProcess

class DeadlineReferenceProcess(ReferenceProcess):
    def call(self,request,timeout=8.):
        if timeout<=0:raise ValueError('Positive execution deadline required')
        result=Queue()
        def read():
            try:result.put((True,ReferenceProcess.call(self,request)))
            except BaseException as exc:result.put((False,exc))
        reader=Thread(target=read,daemon=True);reader.start()
        try:ok,value=result.get(timeout=timeout)
        except Empty:
            process=self.process
            if process is None:raise RuntimeError('Source worker was not established before deadline')
            was_live=process.poll() is None
            if was_live:process.kill()
            code=process.wait(timeout=5.)
            reader.join(timeout=5.)
            if reader.is_alive():raise RuntimeError('Source reader did not finish after worker termination')
            self.process=None
            return {'reference_timeout':{'seconds':timeout,'exit_code':code,
                'worker_pid':process.pid,'intentionally_terminated':was_live}}
        reader.join()
        if ok:return value
        raise value
