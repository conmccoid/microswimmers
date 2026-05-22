# microswimmers
QCAM II working group on microswimmers

## Objectives

- [ ] Distribution of swimmers in FDM
  - [ ] Task 1
  - [ ] Task 2
  - [ ] Task 3
- [ ] Effective diffusivity
- [ ] Rotlets
- [ ] Viscoelastic networks

## Background info

## How-to's

### Installing the code base

### Running the code base

### Exporting MATLAB data to ParaView

1. Delete previous files
```matlab
if isfile('particles.h5')
    delete('particles.h5')
end
```

2. Write raw data to an HDF5
```matlab
for k=1:nstore:Nt
    kk=round((k-1)/nstore)+1;
    coords = sprintf('/coords/step_%04d',kk);
    vel = sprintf('/orientation/step_%04d',kk);
    h5create('particles.h5',coords,[Nfd,3])
    h5write('particles.h5',coords,squeeze(x_store(:,:,kk)))
    h5create('particles.h5',vel,[Nfd,3])
    h5write('particles.h5',vel,squeeze(beta_store(:,:,kk)))
end
```

3. Write a descriptor to an XDMF
```matlab
fileID = fopen('particles.xmf','w');
fprintf(fileID,'<?xml version="1.0" ?>\n');
fprintf(fileID,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fileID,'<Xdmf Version="3.0">\n');
fprintf(fileID,'<Domain>\n');
fprintf(fileID,'    <Grid Name="points" GridType="Collection" CollectionType="Temporal">\n');

for k=1:nstore:Nt
    kk=round((k-1)/nstore)+1;
    coords = sprintf('/coords/step_%04d',kk);
    vel = sprintf('/orientation/step_%04d',kk);
    fprintf(fileID,'        <Grid Name="step_%04d">\n',kk);
    fprintf(fileID,'            <Time Value="%.6f"/>\n', k*dt);
    fprintf(fileID,'            <Topology TopologyType="Polyvertex" NumberOfElements="%d">\n',Nfd);
    fprintf(fileID,'            </Topology>\n');
    fprintf(fileID,'            <Geometry GeometryType="XYZ">\n');
    fprintf(fileID,'                <DataItem Format="HDF" DataType="Float" Dimensions="3 %d">\n',Nfd);
    fprintf(fileID,'                    particles.h5:%s\n', coords);
    fprintf(fileID,'                </DataItem>\n');
    fprintf(fileID,'            </Geometry>\n');
    fprintf(fileID,'            <Attribute Name="Orientation" AttributeType="Vector" Center="Node">\n');
    fprintf(fileID,'                <DataItem Format="HDF" Dimensions="3 %d">\n',Nfd);
    fprintf(fileID,'                    particles.h5:%s\n',vel);
    fprintf(fileID,'                </DataItem>\n');
    fprintf(fileID,'            </Attribute>\n');
    fprintf(fileID,'        </Grid>\n');
end
fprintf(fileID,'    </Grid>\n');
fprintf(fileID,'</Domain>\n');
fprintf(fileID,'</Xdmf>\n');
fclose(fileID);
```

### Visualizing in ParaView

Use the 'Glyph' filter.