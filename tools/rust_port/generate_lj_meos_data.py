"""Extract the two LJ MEOS fits and eight noble-gas parameter sets."""
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
source=(ROOT/'src/multiparameter_lj.f90').read_text()

def array(name):
    match=re.search(r'::\s*'+name+r'\s*=\s*\(/(.*?)/\)',source,re.S|re.I)
    assert match,name
    return [float(v.strip()) for v in match[1].replace('&','').split(',')]

def scalar(block,name):
    match=re.search(r'this%'+name+r'\s*=\s*([\d.eE+\-]+)',block,re.I)
    assert match,name
    return float(match[1])

lines=['//! Generated from src/multiparameter_lj.f90 by generate_lj_meos_data.py.',
       '//! Copyright thermotools and Contributors. Apache-2.0.',
       '#![allow(clippy::excessive_precision)]',
       'use crate::multiparameter_data::{AncillarySatDensity,IdealGasHelmholtz,MultiparameterRecord,PolTerm,ExpTerm,ExpExpTerm};',
       '#[derive(Debug,Clone,Copy)]',
       'pub(crate) struct Gas {pub component: &\'static str,pub truncated:bool,pub epsilon:f64,pub sigma:f64}',
       'pub(crate) const GASES: [Gas;8] = [']
for label,truncated in [('LJ',False),('LJTS',True)]:
    block=re.search(r'::\s*'+label+r'array\s*=\s*\(/(.*?)/\)',source,re.S|re.I).group(1)
    records=re.findall(r'comp_name="([^"]+)",\s*eps_divk=([\d.]+),\s*sigma=([\deE.\-]+)',block)
    assert len(records)==4
    for component,eps,sigma in records:
        lines.append(f'Gas {{ component:"{component}",truncated:{str(truncated).lower()},epsilon:{float(eps)!r},sigma:{float(sigma)!r} }},')
lines.append('];')
for label in ['LJ','LJTS']:
    suffix=label.lower()
    block=re.search(r'subroutine\s+init_'+label+r'\s*\(this\)(.*?)end subroutine\s+init_'+label,source,re.S|re.I).group(1)
    for family,fields,typ in [('pol',['N','t','d'],'PolTerm'),('exp',['N','t','d','l'],'ExpTerm'),('expexp',['N','t','d','eta','beta','gam','eps'],'ExpExpTerm')]:
        arrays=[array(f'{field}_{family}_{suffix}') for field in fields]
        assert len(set(map(len,arrays)))==1
        lines.append(f'static {label}_{family.upper()}: [{typ};{len(arrays[0])}] = [')
        for row in zip(*arrays):
            fields_out=[]
            for field,value in zip(fields,row):
                out_name={'N':'n','gam':'gamma','eps':'epsilon'}.get(field,field)
                fields_out.append(f'{out_name}:{int(value) if field in ["d","l"] else repr(value)}')
            lines.append(typ+'{'+','.join(fields_out)+'},')
        lines.append('];')
    for name in ['N_liqsat','expo_liqsat','N_vapsat','expo_vapsat']:
        data=array(f'{name}_{suffix}')
        lines.append(f'static {label}_{name.upper()}: [f64;{len(data)}]={data!r};')
    # The LJTS initial ideal constant is shifted at runtime in the same
    # sequence as init_LJTS. The source expression is deliberately retained.
    if label=='LJTS':
        expression=re.search(r'::\s*a_ljts\s*=\s*\(/(.*?)/\)',source,re.S|re.I).group(1)
        expressions=[v.strip() for v in expression.replace('&','').split(',')]
        assert len(expressions)==2 and all(re.fullmatch(r'[\d. +*/()eE\-]+',v) for v in expressions)
        # Keep integral float literals in the generator's existing Rust spelling.
        expressions=[v.split('.')[0]+'.' if re.fullmatch(r'[+\-]?\d+\.0+',v) else v for v in expressions]
        a1,a2=expressions
    else:
        a1,a2=map(repr,array('a_lj'))
    fields={name:scalar(block,name) for name in ['tc','pc','rc','acf','t_triple','p_triple','rhoLiq_triple','rhoVap_triple','molarMass','Rgas_fit','maxT','maxP']}
    lines.append(f'pub(crate) const {label}_LIMITS:(f64,f64)=({fields["maxT"]!r},{fields["maxP"]!r});')
    lines.append(f'pub(crate) const {label}_TEMPLATE:MultiparameterRecord=MultiparameterRecord{{ident:"{label}",')
    for field,name in [('tc','tc'),('pc','pc'),('rc','rc'),('acf','acf'),('molar_mass','molarMass'),('t_triple','t_triple'),('p_triple','p_triple'),('rholiq_triple','rhoLiq_triple'),('rhovap_triple','rhoVap_triple'),('rgas_fit','Rgas_fit')]:
        lines.append(f'{field}:{fields[name]!r},')
    lines.extend([f'ideal:IdealGasHelmholtz::Einstein{{c_tau:1.5,a:[{a2},{a1}],b:&[],v:&[]}},',
        f'pol_terms:&{label}_POL,exp_terms:&{label}_EXP,expexp_terms:&{label}_EXPEXP,',
        f'ancillary:AncillarySatDensity::LemmonSpan{{n_liq:&{label}_N_LIQSAT,expo_liq:&{label}_EXPO_LIQSAT,n_vap:&{label}_N_VAPSAT,expo_vap:&{label}_EXPO_VAPSAT}},','};'])
target=ROOT/'rust/src/lj_meos_data.rs'
target.write_text('\n'.join(lines)+'\n')
print(f'Generated {target}: two fits, eight noble-gas parameter sets')
