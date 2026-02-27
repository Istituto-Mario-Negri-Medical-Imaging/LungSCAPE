function [g,nodenums] = binaryImageGraphWeighted(map,conn)
%binaryImageGraphWeighted  Weighted graph of foreground pixel connectivity in 2-D image
%
%   Extension of binaryImageGraph from the Image Graphs toolbox
%   (Steve Eddins, MathWorks File Exchange #53614). This is the 2-D
%   wrapper for binaryImageGraph3Weighted -- it delegates to the 3-D
%   version and strips the z coordinate from the result. Requires the
%   Image Graphs toolbox to be on the MATLAB path (uses requireR2015b).
%
%   g = binaryImageGraphWeighted(MAP,conn)
%
%   [g,nodenums] = binaryImageGraphWeighted(___)
%
%   See also binaryImageGraph, binaryImageGraph3Weighted,
%   plotImageGraph, imageGraph.

%   Based on binaryImageGraph - Copyright 2015 The MathWorks, Inc.
%   Original author: Steve Eddins
%   Weighted extension: Alberto Arrigoni

requireR2015b

% Input argument parsing and validation.
validateattributes(map,{'double','single'},{'2d'});

[g,nodenums] = binaryImageGraph3Weighted(map,conn);
g.Nodes.z = [];
