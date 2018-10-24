function bssePrimedTrackHelp(~)
% Display documentation for Primed Track.
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
%        <Submenu name="Primed Track">
%         <Item name="Documentation" icon="Matlab">
%          <Command>MatlabXT::bssePrimedTrackHelp(%i)</Command>
%         </Item>
%        </Submenu>
%       </Submenu>
%      </Menu>
%    </CustomTools>
%
% Aaron Ponti (BSSE) 2017, 2018

fname = fullfile(fileparts(mfilename('fullpath')), 'bssePrimedTrackHelp.pdf');
system(['start ', fname]);

