"""Line-oriented subprocess for a library that can terminate its caller."""
import argparse
import json
import sys
from reference import load


def main(compiler_directory):
    cubic=load(compiler_directory)
    model,previous=None,None
    for line in sys.stdin:
        request=json.loads(line)
        eos,comps=request['eos'],request['components']
        if previous!=(eos,comps):
            model=cubic(comps,"SRK" if eos=="RK" else eos,
                        alpha="RK" if eos=="RK" else "Classic",volume_shift=False)
            previous=eos,comps
        if request['mode']=='TP':
            result=model.two_phase_tpflash(request['temperature'],request['pressure'],request['feed'])
            values=[int(result.phase),float(result.betaV),float(result.betaL),
                    *map(float,result.x),*map(float,result.y)]
            swappable=False
            if result.phase==model.TWOPH:
                swappable=True
                for composition in [result.x,result.y]:
                    vl=model.specific_volume(request['temperature'],request['pressure'],composition,model.LIQPH)[0]
                    vv=model.specific_volume(request['temperature'],request['pressure'],composition,model.VAPPH)[0]
                    swappable &= abs(vl-vv)<=1e-10*max(abs(vl),abs(vv))
        else:
            raise ValueError('unsupported oracle mode')
        print('THERMOPACK_RESULT='+json.dumps({'values':values,'phase_permutation_allowed':swappable}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    main(parser.parse_args().compiler_directory)
