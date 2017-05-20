using System;

namespace AstroClock
{
    public class AstroClock
    {
        private double latitude;
        private double longitude;
        private string timeZoneId;

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="tzId">TimeZoneId used to convert UTC times to local times</param>
        /// <param name="lat">Latitude in decimal format (like: 47.699732)</param>
        /// <param name="lon">Longitude in decimal format (like: 122.235731)</param>
        public AstroClock(string tzId, double lat, double lon)
        {
            this.timeZoneId = tzId;
            this.latitude = lat;
            this.longitude = lon;
        }

        public AstroClock(string tzId, double degreesLat, double minutesLat, double secondsLat, double degreesLon, double minutesLon, double secondsLon)
        {
            bool latSignChange = false;
            bool lonSignChange = false;

            this.timeZoneId = tzId;

            if (degreesLat < 0.0)
            {
                degreesLat = Math.Abs(degreesLat);
                latSignChange = true;
            }

            if (degreesLon < 0.0)
            {
                degreesLon = Math.Abs(degreesLon);
                lonSignChange = true;
            }

            double lat = degreesLat + minutesLat / 60.0 + secondsLat / 3600.0;
            double lon = degreesLon + minutesLon / 60.0 + secondsLon / 3600.0;

            if (latSignChange)
                this.latitude = -lat;
            else
                this.latitude = lat;

            if (lonSignChange)
                this.longitude = -lon;
            else
                this.longitude = lon;
        }

        /// <summary>
        /// Calculates sunrise time for the date specified in by the input parameter.
        /// </summary>
        /// <param name="dt">DateTime representing the date to calculate sunrise for</param>
        /// <returns>DateTime object containing the calculated sunrise time</returns>
        public DateTime GetSunrise(DateTime dt)
        {
            double jd = JulianDateFromYMD(dt.Year, dt.Month, dt.Day);
            double srUTCMinutes = SunriseUTC(jd);

            DateTime srUTCDateTime = DateTimeFromMinutes(srUTCMinutes, jd);
            TimeZoneInfo tzInfo = TimeZoneInfo.FindSystemTimeZoneById(this.timeZoneId);
            DateTime srLocalTime = TimeZoneInfo.ConvertTimeFromUtc(srUTCDateTime, tzInfo);

            return srLocalTime;
        }

        /// <summary>
        /// Calculates sunset time for the date specified in by the input parameter.
        /// </summary>
        /// <param name="dt">DateTime representing the date to calculate sunset for</param>
        /// <returns>DateTime object containing the calculated sunset time</returns>
        public DateTime GetSunset(DateTime dt)
        {
            double jd = JulianDateFromYMD(dt.Year, dt.Month, dt.Day);
            double ssUTCMinutes = SunsetUTC(jd);

            DateTime ssUTCDateTime = DateTimeFromMinutes(ssUTCMinutes, jd);
            TimeZoneInfo tzInfo = TimeZoneInfo.FindSystemTimeZoneById(this.timeZoneId);
            DateTime ssLocalTime = TimeZoneInfo.ConvertTimeFromUtc(ssUTCDateTime, tzInfo);

            return ssLocalTime;
        }

        /// <summary>
        /// Calculates when the next sunrise event will occur from "now".
        /// </summary>
        /// <returns>DateTime object containg the date and time of the next sunrise event</returns>
        public DateTime GetNextSunrise()
        {
            DateTime dt = DateTime.Now;
            TimeZoneInfo tzInfo = TimeZoneInfo.FindSystemTimeZoneById(this.timeZoneId);
            DateTime srUTCDateTime;
            DateTime srLocalTime;

            double jd = JulianDateFromYMD(dt.Year, dt.Month, dt.Day);
            double srUTCMinutes = SunriseUTC(jd);

            srUTCDateTime = DateTimeFromMinutes(srUTCMinutes, jd);
            srLocalTime = TimeZoneInfo.ConvertTimeFromUtc(srUTCDateTime, tzInfo);

            Int32 c = srLocalTime.CompareTo(DateTime.Now);

            if (srLocalTime.CompareTo(DateTime.Now) < 0)
            {
                jd = JulianDateFromYMD(dt.Year, dt.Month, dt.Day + 1);
                srUTCMinutes = SunriseUTC(jd);
                srUTCDateTime = DateTimeFromMinutes(srUTCMinutes, jd);
                srLocalTime = TimeZoneInfo.ConvertTimeFromUtc(srUTCDateTime, tzInfo);
            }

            return srLocalTime;
        }

        /// <summary>
        /// Calculates when the next sunset event will occur from "now".
        /// </summary>
        /// <returns>DateTime object containg the date and time of the next sunset event</returns>
        public DateTime GetNextSunset()
        {
            DateTime dt = DateTime.Now;
            TimeZoneInfo tzInfo = TimeZoneInfo.FindSystemTimeZoneById(this.timeZoneId);
            DateTime ssUTCDateTime;
            DateTime ssLocalTime;

            double jd = JulianDateFromYMD(dt.Year, dt.Month, dt.Day);
            double ssUTCMinutes = SunsetUTC(jd);

            ssUTCDateTime = DateTimeFromMinutes(ssUTCMinutes, jd);
            ssLocalTime = TimeZoneInfo.ConvertTimeFromUtc(ssUTCDateTime, tzInfo);

            if (ssLocalTime.CompareTo(DateTime.Now) < 0)
            {
                jd = JulianDateFromYMD(dt.Year, dt.Month, dt.Day + 1);
                ssUTCMinutes = SunsetUTC(jd);
                ssUTCDateTime = DateTimeFromMinutes(ssUTCMinutes, jd);
                ssLocalTime = TimeZoneInfo.ConvertTimeFromUtc(ssUTCDateTime, tzInfo);
            }

            return ssLocalTime;
        }

        /// <summary>
        /// Converts radian to degrees.
        /// </summary>
        /// <param name="radian">Radian value to be converted</param>
        /// <returns>Degrees</returns>
        private double RadianToDegrees(double radian) 
        {
            return (180.0 * radian / Math.PI);
        }
     
        /// <summary>
        /// Converts degrees to radian.
        /// </summary>
        /// <param name="degrees">Degree value to convert</param>
        /// <returns>Radian</returns>
        private double DegreesToRadian(double degrees) 
        {
            return (Math.PI * degrees / 180.0);
        }
     
        /// <summary>
        /// Finds numerical day-of-year from month, day and isLeaYear info.
        /// </summary>
        /// <param name="month">January = 1</param>
        /// <param name="day">1 - 31</param>
        /// <param name="isLeapYear">1 if leap year, 0 if not</param>
        /// <returns>The numerical day of year</returns>
        private double DayOfYear(Int32 month, Int32 day, bool isLeapYear) 
        {
            double c = (isLeapYear ? 1.0 : 2.0);

            return Math.Floor((275.0 * month)/9.0) - c * Math.Floor((month + 9.0)/12.0) + day - 30.0;
        }

        /// <summary>
        /// Derives weekday from Julian Day.
        /// </summary>
        /// <param name="jd">Julian day to be mapped to a day of the week</param>
        /// <returns>String value representing the day of the week</returns>
        private string DayOfWeek(double jd)
        {
            double A = (jd + 1.5) % 7;

            return (A==0)?"Sunday":(A==1)?"Monday":(A==2)?"Tuesday":(A==3)?"Wednesday":(A==4)?"Thursday":(A==5)?"Friday":"Saturday";
        }

        /// <summary>
        /// Calculate julian day from calendar day.
        /// </summary>
        /// <param name="year">4 digit year</param>
        /// <param name="month">January = 1</param>
        /// <param name="day">1 - 31</param>
        /// <returns>The Julian day corresponding to the date.  Value represents the start of the day.</returns>
        private double JulianDateFromYMD(Int32 year, Int32 month, Int32 day)
        {
            if (month <= 2)
            {
                year -= 1;
                month += 12;
            }

            double A = Math.Floor(year/100.0);
            double B = 2 - A + Math.Floor(A/4);
            double JD = Math.Floor(365.25*(year + 4716)) + Math.Floor(30.6001*(month+1)) + day + B - 1524.5;
            
            return JD;
        }

        /// <summary>
        /// Calendar day (minus year) from Julian Day.
        /// </summary>
        /// <param name="yearRef">year (out) of the julian date</param>
        /// <param name="monthRef">month (out) of the julian date</param>
        /// <param name="dayRef">day (out) of the julian date</param>
        /// <param name="jd">Julian Day</param>
        private void YMDFromJulianDate(ref Int32 yearRef, ref Int32 monthRef, ref Int32 dayRef, double jd)
        {
            double z = Math.Floor(jd + 0.5);
            double f = (jd + 0.5) - z;
            double A;
     
            if (z < 2299161)
            {
                A = z;
            }
            else
            {
                double alpha = Math.Floor((z - 1867216.25)/36524.25);
                A = z + 1 + alpha - Math.Floor(alpha/4);
            }
     
            double B = A + 1524;
            double C = Math.Floor((B - 122.1)/365.25);
            double D = Math.Floor(365.25 * C);
            double E = Math.Floor((B - D)/30.6001);
     
            double day = B - D - Math.Floor(30.6001 * E) + f;
            double month = (E < 14) ? E - 1 : E - 13;
            double year = (month > 2) ? C - 4716 : C - 4715;

            dayRef = (Int32)day;
            monthRef = (Int32)month;
            yearRef = (Int32)year;
        }

        /// <summary>
        /// Convert Julian Day to centuries since J2000.0.
        /// </summary>
        /// <param name="jd">The Julian Day to convert</param>
        /// <returns>The value corresponding to the Julian Day</returns>
        private double CalcTimeJulianCentury(double jd)
        {
            return (jd - 2451545.0)/36525.0;
        }

        /// <summary>
        /// Convert centuries since J2000.0 to Julian Day.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>The Julian Day corresponding to the t value</returns>
        private double JulianDateFromJulianCentury(double t)
        {
            return (t * 36525.0 + 2451545.0);
        }

        /// <summary>
        /// Calculate the Geometric Mean Longitude of the Sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>The Geometric Mean Longitude of the Sun in degrees</returns>
        private double GeoMeanLongSun(double t)
        {
            double L0 = 280.46646 + t * (36000.76983 + 0.0003032 * t);
            
            while (L0 > 360.0)
            {
                L0 -= 360.0;
            }
            
            while (L0 < 0.0)
            {
                L0 += 360.0;
            }
            
            return L0;
        }

        /// <summary>
        /// Calculate the Geometric Mean Anomaly of the Sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>The Geometric Mean Anomaly of the Sun in degrees</returns>
        private double GeoMeanAnomalySun(double t)
        {
            return (357.52911 + t * (35999.05029 - 0.0001537 * t));
        }

        /// <summary>
        /// Calculate the eccentricity of earth's orbit.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>The unitless eccentricity</returns>
        private double EccentricityEarthOrbit(double t)
        {
            return (0.016708634 - t * (0.000042037 + 0.0000001267 * t));
        }
     
        /// <summary>
        /// Calculate the equation of center for the sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Center in degrees</returns>
        private double SunEqOfCenter(double t)
        {
            double m = GeoMeanAnomalySun(t);
            double mrad = DegreesToRadian(m);
            double sinm = Math.Sin(mrad);
            double sin2m = Math.Sin(mrad+mrad);
            double sin3m = Math.Sin(mrad+mrad+mrad);
            double C = sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289;
            
            return C;
        }

        /// <summary>
        /// Calculate the true longitude of the sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Sun's true longitude in degrees</returns>
        private double SunTrueLong(double t)
        {
            return (GeoMeanLongSun(t) + SunEqOfCenter(t));
        }

        /// <summary>
        /// Calculate the true anamoly of the sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Sun's true anamoly in degrees</returns>
        private double SunTrueAnomaly(double t)
        {
            return (GeoMeanAnomalySun(t) + SunEqOfCenter(t));
        }
     
        /// <summary>
        /// Calculate the distance to the sun in Astronomical Units (AU).
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Sun radius vector in AUs</returns>
        private double SunRadVector(double t)
        {
            double v = SunTrueAnomaly(t);
            double e = EccentricityEarthOrbit(t);
            double R = (1.000001018 * (1 - e * e)) / (1 + e * Math.Cos(DegreesToRadian(v)));
            
            return R;
        }

        /// <summary>
        /// Calculate the apparent longitude of the sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Sun's apparent longitude in degrees</returns>
        private double SunApparentLong(double t)
        {
            double o = SunTrueLong(t);
            double omega = 125.04 - 1934.136 * t;
            double lambda = o - 0.00569 - 0.00478 * Math.Sin(DegreesToRadian(omega));
            
            return lambda;
        }

        /// <summary>
        /// Calculate the mean obliquity of the ecliptic.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Mean obliquity in degrees</returns>
        private double MeanObliquityOfEcliptic(double t)
        {
            double seconds = 21.448 - t*(46.8150 + t*(0.00059 - t*(0.001813)));
            double e0 = 23.0 + (26.0 + (seconds/60.0))/60.0;
            
            return e0;
        }

        /// <summary>
        /// Calculate the corrected obliquity of the ecliptic.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Corrected obliquity in degrees</returns>
        private double ObliquityCorrection(double t)
        {
            double e0 = MeanObliquityOfEcliptic(t);
            double omega = 125.04 - 1934.136 * t;
            double e = e0 + 0.00256 * Math.Cos(DegreesToRadian(omega));
            
            return e;
        }

        /// <summary>
        /// Calculate the right ascension of the sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Sun's right ascension in degrees</returns>
        private double SunRightAscension(double t)
        {
            double e = ObliquityCorrection(t);
            double lambda = SunApparentLong(t);
            double tananum = (Math.Cos(DegreesToRadian(e)) * Math.Sin(DegreesToRadian(lambda)));
            double tanadenom = (Math.Cos(DegreesToRadian(lambda)));
            double alpha = RadianToDegrees(Math.Atan2(tananum, tanadenom));
            
            return alpha;
        }

        /// <summary>
        /// Calculate the declination of the sun.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Sun's declination in degrees</returns>
        private double SunDeclination(double t)
        {
            double e = ObliquityCorrection(t);
            double lambda = SunApparentLong(t);
            double sint = Math.Sin(DegreesToRadian(e)) * Math.Sin(DegreesToRadian(lambda));
            double theta = RadianToDegrees(Math.Asin(sint));
            
            return theta;
        }

        /// <summary>
        /// Calculate the difference between true solar time and mean solar time.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Equation of time in minutes of time</returns>
        private double EquationOfTime(double t)
        {
            double epsilon = ObliquityCorrection(t);
            double l0 = GeoMeanLongSun(t);
            double e = EccentricityEarthOrbit(t);
            double m = GeoMeanAnomalySun(t);
            double y = Math.Tan(DegreesToRadian(epsilon)/2.0);
            
            y *= y;
     
            double sin2l0 = Math.Sin(2.0 * DegreesToRadian(l0));
            double sinm   = Math.Sin(DegreesToRadian(m));
            double cos2l0 = Math.Cos(2.0 * DegreesToRadian(l0));
            double sin4l0 = Math.Sin(4.0 * DegreesToRadian(l0));
            double sin2m  = Math.Sin(2.0 * DegreesToRadian(m));
            double Etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;
     
            return (RadianToDegrees(Etime)*4.0);
        }

        /// <summary>
        /// Calculate the hour angle of the sun at sunrise for the latitude.
        /// </summary>
        /// <param name="solarDec">Declination angle of sun in degrees</param>
        /// <returns>Hour angle of sunrise in radians</returns>
        private double HourAngleSunrise(double solarDec)
        {
            double latRad = DegreesToRadian(this.latitude);
            double sdRad  = DegreesToRadian(solarDec);
            double HA = (Math.Acos(Math.Cos(DegreesToRadian(90.833))/(Math.Cos(latRad)*Math.Cos(sdRad))-Math.Tan(latRad) * Math.Tan(sdRad)));
     
            return HA;
        }

        /// <summary>
        /// Calculate the hour angle of the sun at sunset for the latitude.
        /// </summary>
        /// <param name="solarDec">Declination angle of sun in degrees</param>
        /// <returns>Hour angle of sunset in radians</returns>
        private double HourAngleSunset(double solarDec)
        {
            double latRad = DegreesToRadian(this.latitude);
            double sdRad = DegreesToRadian(solarDec);
            double HA = (Math.Acos(Math.Cos(DegreesToRadian(90.833)) / (Math.Cos(latRad) * Math.Cos(sdRad)) - Math.Tan(latRad) * Math.Tan(sdRad)));

            return -HA;
        }

        /// <summary>
        /// Calculate the Universal Coordinated Time (UTC) of sunrise for the 
        /// given day at the given location on earth.
        /// </summary>
        /// <param name="jd">Julian day</param>
        /// <returns>Time in minutes from zero Z</returns>
        private double SunriseUTC(double jd)
        {
            double t = CalcTimeJulianCentury(jd);
     
            // Find the time of solar noon at the location, and use that
            // declination. This is better than start of the Julian day
     
            double noonmin = SolarNoonUTC(t);
            double tnoon = CalcTimeJulianCentury(jd+noonmin/1440.0);
     
            // First pass to approximate sunrise (using solar noon)
     
            double eqTime = EquationOfTime(tnoon);
            double solarDec = SunDeclination(tnoon);
            double hourAngle = HourAngleSunrise(solarDec);
     
            double delta = longitude - RadianToDegrees(hourAngle);
            double timeDiff = 4 * delta;               // in minutes of time
            double timeUTC = 720 + timeDiff - eqTime;  // in minutes
     
            // Second pass includes fractional jday in gamma calc
     
            double newt = CalcTimeJulianCentury(JulianDateFromJulianCentury(t) + timeUTC/1440.0); 
            
            eqTime = EquationOfTime(newt);
            solarDec = SunDeclination(newt);
            hourAngle = HourAngleSunrise(solarDec);
            delta = longitude - RadianToDegrees(hourAngle);
            timeDiff = 4 * delta;
            timeUTC = 720 + timeDiff - eqTime; // in minutes
     
            return timeUTC;
        }

        /// <summary>
        /// Calculate the Universal Coordinated Time (UTC) of solar noon for 
        /// the given day at the given location on earth.
        /// </summary>
        /// <param name="t">Number of Julian centuries since J2000.0</param>
        /// <returns>Time in minutes from zero Z</returns>
        private double SolarNoonUTC(double t)
        {
            // First pass uses approximate solar noon to calculate eqtime
            double tnoon = CalcTimeJulianCentury(JulianDateFromJulianCentury(t) + this.longitude/360.0);
            double eqTime = EquationOfTime(tnoon);
            double solarNoonUTC = 720 + (this.longitude * 4) - eqTime; // min
     
            double newt = CalcTimeJulianCentury(JulianDateFromJulianCentury(t) -0.5 + solarNoonUTC/1440.0); 
     
            eqTime = EquationOfTime(newt);
            // double solarNoonDec = SunDeclination(newt);
            solarNoonUTC = 720 + (longitude * 4) - eqTime; // min
            
            return solarNoonUTC;
        }

        /// <summary>
        /// Calculate the Universal Coordinated Time (UTC) of sunset for the 
        /// given day at the given location on earth.
        /// </summary>
        /// <param name="jd">Julian day</param>
        /// <returns>Time in minutes from zero Z</returns>
        private double SunsetUTC(double jd)
        {
            double t = CalcTimeJulianCentury(jd);
     
            // *** Find the time of solar noon at the location, and use
            //     that declination. This is better than start of the 
            //     Julian day
     
            double noonmin = SolarNoonUTC(t);
            double tnoon = CalcTimeJulianCentury (jd+noonmin/1440.0);
     
            // First calculates sunrise and approx length of day
     
            double eqTime = EquationOfTime(tnoon);
            double solarDec = SunDeclination(tnoon);
            double hourAngle = HourAngleSunset(solarDec);
     
            double delta = longitude - RadianToDegrees(hourAngle);
            double timeDiff = 4 * delta;
            double timeUTC = 720 + timeDiff - eqTime;
     
            // first pass used to include fractional day in gamma calc
     
            double newt = CalcTimeJulianCentury(JulianDateFromJulianCentury(t) + timeUTC/1440.0); 
            
            eqTime = EquationOfTime(newt);
            solarDec = SunDeclination(newt);
            hourAngle = HourAngleSunset(solarDec);
     
            delta = longitude - RadianToDegrees(hourAngle);
            timeDiff = 4 * delta;
            timeUTC = 720 + timeDiff - eqTime; // in minutes
     
            return timeUTC;
        }

        /// <summary>
        /// Convert time of day in minutes and the julian date to a .Net DateTime
        /// object.  The time is in UTC and will need to be converted to local
        /// time outside of this routine if that is needed.
        /// </summary>
        /// <param name="minutes">Time of day in minutes</param>
        /// <param name="jd">Julian day</param>
        /// <returns>DateTime from minues (for the time) and jd (for the date)</returns>
        private DateTime DateTimeFromMinutes(double minutes, double jd)
        {
            double floatHour = minutes / 60.0;
            double hour = Math.Floor(floatHour);
            double floatMinute = 60.0 * (floatHour - Math.Floor(floatHour));
            double minute = Math.Floor(floatMinute);
            double floatSec = 60.0 * (floatMinute - Math.Floor(floatMinute));
            double second = Math.Floor(floatSec + 0.5);
     
            minute += (second >= 30)? 1 : 0;
     
            if (minute >= 60) 
            {
                minute -= 60;
                hour ++;
            }
     
            if (hour > 23) 
            {
                hour -= 24;
                jd += 1.0;
            }
     
            if (hour < 0)
            {
                hour += 24;
                jd -= 1.0;
            }

            Int32 year = 0;
            Int32 month = 0;
            Int32 day = 0;

            YMDFromJulianDate(ref year, ref month, ref day, jd);

            return new DateTime(year, month, day, (Int32)hour, (Int32)minute, 0);
        }
    }
}
