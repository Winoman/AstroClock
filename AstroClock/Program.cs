using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AstroClock
{
    class Program
    {
        static void Main(string[] args)
        {
            //AstroClock ac = new AstroClock("India Standard Time", 12.0, 58.0, 0.0, -77.0, 34.0, 0.0);
            AstroClock ac = new AstroClock("Pacific Standard Time", 47.699732, 122.235731);
            DateTime dt = new DateTime(2017, 5, 20, 3, 15, 0);

            DateTime dtSunrise = ac.GetSunrise(DateTime.Now);
            DateTime dtSunset = ac.GetSunset(DateTime.Now);
            DateTime dtNextSunrise = ac.GetNextSunrise();
            DateTime dtNextSunset = ac.GetNextSunset();

            Console.WriteLine("Sunrise: {0}", dtSunrise.ToString("g"));
            Console.WriteLine("Sunset:  {0}", dtSunset.ToString("g"));
            Console.WriteLine("Next Sunrise: {0}", dtNextSunrise.ToString("g"));
            Console.WriteLine("Next Sunset:  {0}", dtNextSunset.ToString("g"));

            DateTime dtSunriseJune = ac.GetSunrise(dt);
            DateTime dtSunsetJune = ac.GetSunset(dt);

            Console.WriteLine("Sunrise: {0}", dtSunriseJune.ToString("g"));
            Console.WriteLine("Sunset:  {0}", dtSunsetJune.ToString("g"));
            
            Console.ReadLine();
        }
    }
}
