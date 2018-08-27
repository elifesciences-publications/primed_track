function bsseLineageTracerHelp(~)
% Display documentation for the Lineage Tracer.
%
% Code for the paper:
%
% Welling et al. "High fidelity lineage tracing in mouse pre-implantation 
% embryos using primed conversion of photoconvertible proteins".
%
% This Imaris XTension required IceImarisConnector to run.
% See: https://github.com/aarpon/IceImarisConnector
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="BSSE">
%        <Submenu name="Lineage Tracer">
%         <Item name="Documentation" icon="Matlab">
%          <Command>MatlabXT::bsseLineageTracerHelp(%i)</Command>
%         </Item>
%        </Submenu>
%       </Submenu>
%      </Menu>
%    </CustomTools>
%
% Aaron Ponti (BSSE) 2017, 2018

fname = fullfile(fileparts(mfilename('fullpath')), 'bsseLineageTracerHelp.pdf');
system(['start ', fname]);

