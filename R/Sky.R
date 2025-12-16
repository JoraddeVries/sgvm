# ==========================
# Astronomical calculations
# ==========================

# Declination angle (radians)
declination <- function(DOY) {
  23.45 * pi / 180 * sin(2 * pi / 365 * (DOY - 81))
}

# Solar angles: returns list(cos_theta, theta, phi)
solar_angles <- function(lat, t, dec) {
  h <- (t - 12) * 15 * pi / 180      # hour angle in radians
  cos_theta <- cos(lat) * cos(dec) * cos(h) + sin(lat) * sin(dec)
  theta <- acos(cos_theta)            # zenith angle
  cos_phi <- (sin(dec)*cos(lat) - cos(h)*cos(dec)*sin(lat)) / sin(theta)
  cos_phi <- pmax(pmin(cos_phi, 1), -1)
  phi <- acos(cos_phi)
  phi <- ifelse(h < 0, phi, 2*pi - phi)
  list(cos_theta = cos_theta, theta = theta, phi = phi)
}

# Extraterrestrial radiation (W/m^2)
extraterrestrial <- function(DOY) {
  1367.7 * (1 + 0.033 * cos(DOY / 365 * 2 * pi))
}

# Air mass
air_mass <- function(cos_theta, theta) {
  1 / (cos_theta + 0.50572 * (96.07995 - theta*180/pi)^(-1.6354))
}

# Day length (hours)
day_length <- function(lat, dec) {
  cosSunset <- -tan(lat) * tan(dec)
  sunset <- ifelse(cosSunset > 1, 0,
                   ifelse(cosSunset < -1, pi, acos(cosSunset)))
  2 * sunset * 180 / pi / 15
}

# Time of day (hours)
timeOfDay <- function(f, DL) {
  tsr <- 12 - DL / 2
  tsr + f * DL
}

# ==========================
# Clear-sky radiation
# ==========================

clear_sky <- function(lat, DOY, f, altitude = 0, TL = 4) {
  stopifnot(abs(lat) <= pi/2, 0 <= f, f <= 1, DOY > 0, DOY <= 365)
  Io <- extraterrestrial(DOY)
  dec <- declination(DOY)
  DL <- day_length(lat, dec)
  t <- timeOfDay(f, DL)
  angles <- solar_angles(lat, t, dec)
  cos_theta <- angles$cos_theta
  theta <- angles$theta
  phi <- angles$phi
  am <- air_mass(cos_theta, theta)
  
  fh1 <- exp(-altitude / 8000)
  fh2 <- exp(-altitude / 1250)
  a1 <- 5.09e-5 * altitude + 0.868
  a2 <- 3.92e-5 * altitude + 0.0387
  Ig <- Io * cos_theta * a1 * exp(-a2 * am * (fh1 + fh2 * (TL - 1)))
  b <- 0.664 + 0.163 / fh1
  Idir <- Io * cos_theta * b * exp(-0.09 * am * (TL - 1))
  Idif <- Ig - Idir
  list(Ig = Ig, Idir = Idir, Idif = Idif, theta = theta, phi = phi, dayLength = DL)
}

# ==========================
# Daily radiation (Spitters 1986)
# ==========================

daily_radiation <- function(lat, DOY, Igd = NULL) {
  stopifnot(abs(lat) <= pi/2, DOY > 0, DOY <= 365)
  dec <- declination(DOY)
  DL <- day_length(lat, dec)
  int <- 3600 * (DL*sin(lat)*sin(dec) + (24/pi)*cos(lat)*cos(dec)*sqrt(1 - (tan(lat)^2)*(tan(dec)^2)))
  Iod <- extraterrestrial(DOY) * int
  if (is.null(Igd)) {
    Igd <- 0.75 * Iod
    Idif <- 0.23 * Igd
  } else {
    tau <- Igd/Iod
    fdif <- if (tau >= 0.75) 0.23 else if (tau > 0.35) 1.33 - 1.46*tau else if (tau > 0.07) 1 - 2.3*(tau - 0.07)^2 else 1
    Idif <- fdif * Igd
  }
  Idir <- Igd - Idif
  list(Iod = Iod, Igd = Igd, Idif = Idif, Idir = Idir)
}

# ==========================
# Cloudy sky radiation
# ==========================

cloudy_sky <- function(Ig = NULL, Iday = NULL, lat, DOY, f) {
  stopifnot(abs(lat) <= pi/2, 0 <= f, f <= 1, DOY > 0, DOY <= 365)
  stopifnot(!is.null(Ig) || !is.null(Iday))
  Io <- extraterrestrial(DOY)
  dec <- declination(DOY)
  DL <- day_length(lat, dec)
  t <- timeOfDay(f, DL)
  angles <- solar_angles(lat, t, dec)
  theta <- angles$theta
  phi <- angles$phi
  beta <- pi/2 - theta
  
  if (is.null(Ig) && !is.null(Iday)) {
    ft <- Io*cos(theta)/Iday$Iod
    Ig <- ft * Iday$Igd
    Idif <- ft * Iday$Idif
  } else if (!is.null(Ig)) {
    IgIo <- ifelse(beta == 0, 0, min(1, Ig/Io/sin(beta)))
    R <- 0.847 - 1.61*sin(beta) + 1.04*sin(beta)^2
    K <- (1.47 - R)/1.66
    IdfIg <- if (IgIo <= 0.22) 1 else if (IgIo <= 0.35) 1 - 6.4*(IgIo - 0.22)^2 else if (IgIo <= K) 1.47 - 1.66*IgIo else R
    Idif <- Ig * IdfIg
  } else stop("Assign a value to either Ig or Iday, not both.")
  
  list(Ig = Ig, Idir = Ig - Idif, Idif = Idif, theta = theta, phi = phi, dayLength = DL)
}

# ==========================
# Waveband conversion (Bird model)
# ==========================

wavebands_dict <- list(
  direct = list(
    power = list(UV=0.03, blue=0.12, green=0.16, red=0.13, PAR=0.40, NIR=0.56),
    flux  = list(UV=0.03/0.327, blue=0.12/0.263, green=0.16/0.217, red=0.13/0.184, PAR=0.40/0.216, NIR=0.56/0.105)
  ),
  diffuse = list(
    power = list(UV=0.13, blue=0.25, green=0.22, red=0.13, PAR=0.60, NIR=0.27),
    flux  = list(UV=0.13/0.332, blue=0.25/0.265, green=0.22/0.219, red=0.13/0.185, PAR=0.60/0.226, NIR=0.27/0.122)
  )
)

waveband_conversion <- function(Itype = "direct", waveband = "PAR", mode = "power") {
  stopifnot(Itype %in% c("direct","diffuse"), waveband %in% c("PAR","UV","blue","red","green","NIR"),
            mode %in% c("power","flux"))
  wavebands_dict[[Itype]][[mode]][[waveband]]
}
