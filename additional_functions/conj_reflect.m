function cnj = conj_reflect(array)

% Jesse Clark, University College London, 2012

F=ifftshift(fftn(fftshift(array)));

cnj=fftshift(ifftn(ifftshift(conj(F))));


end
