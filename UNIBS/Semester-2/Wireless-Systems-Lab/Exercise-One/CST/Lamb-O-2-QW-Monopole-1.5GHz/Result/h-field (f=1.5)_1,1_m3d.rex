<?xml version="1.0" encoding="UTF-8"?>
<MetaResultFile version="20200629" creator="Solver HFTD - Field 3DFD Monitor">
  <MetaGeometryFile filename="model.gex" lod="1"/>
  <SimulationProperties fieldname="h-field (f=1.5) [1]" fieldtype="H-Field" frequency="1.5" encoded_unit="&amp;U:A^1.:m^-1" fieldscaling="PEAK" dB_Amplitude="20"/>
  <ResultDataType vector="1" complex="1" timedomain="0"/>
  <SimulationDomain min="-90.8140564 -90.8140564 -37.4740562" max="90.8140564 90.8140564 83.0353088"/>
  <PlotSettings Plot="4" ignore_symmetry="0" deformation="0" enforce_culling="0" default_arrow_type="ARROWS"/>
  <Source type="SOLVER"/>
  <SpecialMaterials>
    <Background type="NORMAL"/>
    <Material name="PEC" type="FIELDFREE"/>
  </SpecialMaterials>
  <AuxGeometryFile/>
  <AuxResultFile/>
  <FieldFreeNodes/>
  <SurfaceFieldCoefficients/>
  <UnitCell/>
  <SubVolume/>
  <Units/>
  <ResultGroups num_steps="1" transformation="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1">
    <SharedDataWith treepath="2D/3D Results\Surface Current\surface current (f=1.5) [1]" rexname="h-field (f=1.5)_1,1_m3d_sct.rex"/>
    <Frame index="0">
      <FieldResultFile filename="h-field (f=1.5)_1,1.m3d" type="m3d"/>
    </Frame>
  </ResultGroups>
</MetaResultFile>
