ncfrog <- c(104, 110, 117, 121, 127, 130, 136, 145, 152, 178, 184, 192, 200, 316)
ncfrog.rad <- ncfrog/180*pi

if (require(circular)) {
    ncfrog <- circular(ncfrog, units="degrees", template="geographics")
    ncfrog.rad <- circular(ncfrog.rad, template="geographics")
}
