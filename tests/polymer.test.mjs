import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {buildPolymer, getTemplate, toPDB, toSDF} from '../src/polymer.js';

const templates = JSON.parse(readFileSync(new URL('../data/templates.json', import.meta.url)));
const edge = (a,b) => [a,b].sort((x,y)=>x-y).join(':');

function parsePDB(text) {
  const atoms=[], edges=new Set(), directedCounts=new Map();
  for (const line of text.split('\n')) {
    if (line.startsWith('HETATM')) {
      const q=line.slice(78,80).trim();
      atoms.push({id:Number(line.slice(6,11)), name:line.slice(12,16).trim(), residue:Number(line.slice(22,26)),
        residueName:line.slice(17,20).trim(), element:line.slice(76,78).trim(),
        charge:q?Number(q[0])*(q[1]==='+'?1:-1):0,
        x:Number(line.slice(30,38)), y:Number(line.slice(38,46)), z:Number(line.slice(46,54))});
    }
    if (line.startsWith('CONECT')) {
      const ids=line.slice(6).match(/.{1,5}/g).map(Number);
      assert.ok(ids.length<=5,'CONECT has at most four target fields');
      for (const id of ids.slice(1)) {
        edges.add(edge(ids[0]-1,id-1));
        const key=`${ids[0]-1}:${id-1}`;
        directedCounts.set(key,(directedCounts.get(key)||0)+1);
      }
    }
  }
  const bondOrders=new Map([...edges].map(key=>{
    const [a,b]=key.split(':').map(Number);
    assert.equal(directedCounts.get(`${a}:${b}`),directedCounts.get(`${b}:${a}`));
    return [key,directedCounts.get(`${a}:${b}`)];
  }));
  return {atoms,edges,bondOrders};
}

function parseV3000(text) {
  assert.equal(text.split('\n')[1].slice(20,22),'3D');
  const atoms=[],bonds=[]; let section=''; let counts;
  for (const line of text.split('\n')) {
    if (!line.startsWith('M  V30 ')) continue;
    const fields=line.slice(7).split(/\s+/);
    if (fields[0]==='COUNTS') counts=fields.slice(1,3).map(Number);
    else if (fields[0]==='BEGIN') section=fields[1];
    else if (fields[0]==='END') section='';
    else if (section==='ATOM') atoms.push({id:Number(fields[0]),element:fields[1],
      x:Number(fields[2]),y:Number(fields[3]),z:Number(fields[4]),
      charge:Number(fields.find(x=>x.startsWith('CHG='))?.slice(4)||0)});
    else if (section==='BOND') bonds.push({a:Number(fields[2])-1,b:Number(fields[3])-1,order:Number(fields[1])});
  }
  assert.deepEqual(counts,[atoms.length,bonds.length]);
  return {atoms,bonds};
}

test('all 16 source pairs match recorded SHA256 provenance',()=>{
  assert.equal(Object.keys(templates).length,16);
  for (const template of Object.values(templates)) for (const info of Object.values(template.source)) {
    const bytes=readFileSync(new URL('../'+info.path,import.meta.url));
    assert.equal(createHash('sha256').update(bytes).digest('hex'),info.sha256);
  }
});

for (const template of Object.values(templates)) for (const repeats of [1,2,12,100]) {
  test(`${template.name} DP${repeats}: PDB and V3000 graph/charge/coordinate roundtrips`,()=>{
    const model=buildPolymer(templates,{group:template.group,carbons:template.carbons,repeats});
    const pdb=parsePDB(toPDB(model)),sdf=parseV3000(toSDF(model));
    assert.equal(model.atoms.length,template.atoms.length*repeats-2*repeats+8);
    assert.equal(model.bonds.length,template.bonds.length*repeats-repeats+7);
    assert.equal(model.metadata.removedHydrogens,2*repeats);
    assert.equal(model.metadata.addedCapAtoms,8);
    assert.equal(pdb.atoms.length,model.atoms.length);
    assert.deepEqual(pdb.atoms.map(a=>a.residue),pdb.atoms.map(a=>a.residue).sort((a,b)=>a-b));
    assert.deepEqual(pdb.edges,new Set(model.bonds.map(b=>edge(b.a,b.b))));
    assert.deepEqual(pdb.bondOrders,new Map(model.bonds.map(b=>[edge(b.a,b.b),b.order])));
    assert.deepEqual(sdf.bonds,model.bonds);
    assert.equal(sdf.atoms.reduce((sum,a)=>sum+a.charge,0),repeats);
    for (let i=0;i<model.atoms.length;i++) {
      const original=model.atoms[i],p=pdb.atoms[i],s=sdf.atoms[i];
      for (const key of ['id','name','element','residue','residueName','charge']) assert.equal(p[key],original[key]);
      for (const key of ['id','element','charge']) assert.equal(s[key],original[key]);
      for (const key of ['x','y','z']) {assert.ok(Math.abs(p[key]-original[key])<=0.000501);assert.ok(Math.abs(s[key]-original[key])<=0.000000501);}
    }
  });
}

test('builds are deterministic and do not mutate the templates',()=>{
  const before=JSON.stringify(templates);
  const options={group:'GD',carbons:10,repeats:4};
  assert.deepEqual(buildPolymer(templates,options),buildPolymer(templates,options));
  assert.equal(JSON.stringify(templates),before);
});

test('both main-chain endpoints have new neutral methyl caps, including DP1',()=>{
  for (const template of Object.values(templates)) for (const repeats of [1,2,12,100]) {
    const model=buildPolymer(templates,{group:template.group,carbons:template.carbons,repeats});
    const neighbors=model.atoms.map(()=>[]);
    for (const {a,b} of model.bonds) {neighbors[a].push(b);neighbors[b].push(a);}
    assert.equal(model.atoms.filter(a=>a.role==='endCap').length,8);
    for (let r=1;r<=repeats;r++) {
      const named=Object.fromEntries(model.atoms.map((a,i)=>[a.residue===r?a.name:'',i]));
      assert.equal(neighbors[named.C1].filter(i=>model.atoms[i].element==='H').length,2);
      assert.equal(neighbors[named.C4].filter(i=>model.atoms[i].element==='H').length,0);
      assert.equal(neighbors[named.C5].filter(i=>model.atoms[i].element==='H').length,3);
      assert.equal(model.atoms.some(a=>a.residue===r&&(a.name==='H1'||a.name==='H7')),false);
    }
    for (const [end,number,residue,attachment] of [['head',1,1,model.metadata.units[0].head],['tail',2,repeats,model.metadata.units.at(-1).tail]]) {
      const cap=model.metadata.endCaps[end], c=model.atoms[cap.carbon];
      assert.equal(cap.attachment,attachment);
      assert.equal(c.name,`CM${number}`);
      assert.equal(c.residue,residue);
      assert.equal(c.element,'C');
      assert.equal(c.charge,0);
      assert.deepEqual(new Set(neighbors[cap.carbon]),new Set([attachment,...cap.hydrogens]));
      for (const [i,h] of cap.hydrogens.entries()) {
        assert.equal(model.atoms[h].name,`HM${number}${i+1}`);
        assert.equal(model.atoms[h].element,'H');
        assert.equal(model.atoms[h].charge,0);
        assert.equal(model.atoms[h].endCap,end);
        assert.deepEqual(neighbors[h],[cap.carbon]);
      }
    }
    const elements=atoms=>atoms.reduce((m,a)=>(m[a.element]=(m[a.element]||0)+1,m),{});
    const actual=elements(model.atoms),source=elements(template.atoms);
    assert.equal(actual.C,repeats*source.C+2);
    assert.equal(actual.H,repeats*source.H-2*repeats+6);
    assert.match(toPDB(model),/BOTH MAIN-CHAIN ENDS CAPPED BY NEUTRAL CH3/);
    assert.match(toSDF(model),/>  <END_CAPS>\nNeutral CH3 at both main-chain endpoints/);
  }
});

test('invalid selections and non-finite export coordinates are rejected',()=>{
  for (const repeats of [0,-1,101,1.5,NaN,'2']) assert.throws(()=>buildPolymer(templates,{group:'PA',carbons:4,repeats}),RangeError);
  assert.throws(()=>getTemplate(templates,'NA',4),RangeError);
  assert.throws(()=>getTemplate(templates,'PA',5),RangeError);
  const model=buildPolymer(templates,{group:'PA',carbons:4,repeats:1});
  model.atoms[0].x=NaN;
  assert.throws(()=>toPDB(model),/finite/);
  assert.throws(()=>toSDF(model),/finite/);
});
