"""Private batch worker. Reads its request only after parent containment."""
import json
import re
import sys
from time import monotonic
import traceback

from .atomic_io import write_json
from .data_models import read_json


class Progress:
    def __init__(self,path):
        self.path=path
        self.current={'stage':'importando_workflow','reference_id':None}
        write_json(path,self.current)

    def line(self,line):
        match=re.search(r'\[si-Fi\] (.+?):',line)
        if match:
            self.current.update(reference_id=match[1],stage='si-Fi21')
        elif '[rnai-designer] ' in line:
            self.current['stage']=line.split('[rnai-designer] ',1)[1].strip()
        else:
            return
        try:
            write_json(self.path,self.current)
        except OSError:
            pass  # Diagnostics must not change the scientific workflow's outcome.


class ProgressStream:
    def __init__(self,stream,progress):
        self.stream,self.progress,self.pending=stream,progress,''

    def write(self,text):
        result=self.stream.write(text)
        self.pending+=text
        while '\n' in self.pending:
            line,self.pending=self.pending.split('\n',1)
            self.progress.line(line)
        self.pending=self.pending[-8192:]
        return result

    def __getattr__(self,name):
        return getattr(self.stream,name)


def main(workflow=None):
    request=json.loads(sys.stdin.readline())
    if request.get('request_file'):
        request={**read_json(request['request_file']), 'started':request['started'],
                 'deadline':request.get('deadline')}
    progress=Progress(request['progress'])
    stdout,stderr=sys.stdout,sys.stderr
    sys.stdout,sys.stderr=ProgressStream(stdout,progress),ProgressStream(stderr,progress)
    code=2
    error=None
    try:
        if workflow is None:
            from .__main__ import workflow_command
            workflow=workflow_command
        progress.current['stage']='workflow'
        write_json(progress.path,progress.current)
        code=workflow(request['argv'],input_metadata=request['metadata'])
    except SystemExit as exc:
        code=exc.code if isinstance(exc.code,int) else 2
        error=str(exc)
    except Exception as exc:
        error=str(exc)
        traceback.print_exc()
    finally:
        sys.stdout,sys.stderr=stdout,stderr
        deadline=request.get('deadline')
        write_json(request['receipt'],{'return_code':code,'error':error,
                   'elapsed_s':monotonic()-request['started'],
                   'within_deadline':deadline is None or monotonic()<=deadline})
    return code


if __name__ == '__main__':
    sys.exit(main())
