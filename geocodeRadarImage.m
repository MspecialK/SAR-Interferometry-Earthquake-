function geocodedImage = geocodeRadarImage(sarImage,geocodingInfos)
%GEOCODERADARIMAGE Geocode the SAR image sarImage using the geocoding
%informations inside geocodingInfos (GG)
%

% get image sizes
[Nxsar, Nysar] = size(sarImage);
if Nxsar ~= geocodingInfos.Nsar(1) || Nysar ~= geocodingInfos.Nsar(2)
    warning('Image size does not match');
    return
end

% geocode
geocodedImage = NaN([geocodingInfos.Ngec(1), geocodingInfos.Ngec(2)],'single');
geocodedImage(geocodingInfos.mask) = sarImage(geocodingInfos.indexes);

end