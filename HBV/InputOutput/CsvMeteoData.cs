using System.Globalization;

namespace HBV.InputOutput
{
    public class CsvMeteoData
    {
        public List<DateTime> allDateTimes = new List<DateTime>();
        public List<float> Tmean = new List<float>();
        public List<float> Tdew = new List<float>();
        public List<float> Tground = new List<float>();
        public List<float> RH = new List<float>();
        public List<float> AP = new List<float>();
        public List<float> u = new List<float>();
        public List<float> Rad = new List<float>();
        public List<float> Rad_MJm2d = new List<float>();
        public List<float> Sunshine = new List<float>();
        public List<float> PET = new List<float>();
        public List<float> Precipitation_mm_even = new List<float>();
        public List<float> Discharge = new List<float>();

        public CsvMeteoData(string path)
        {
            using (var reader = new StreamReader(path))
            {
                string headerLine = reader.ReadLine(); // Read and discard the header
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    string[] values = line.Split(',');

                    try
                    {
                        allDateTimes.Add(DateTime.ParseExact(values[0], "MM/dd/yyyy HH:mm", CultureInfo.InvariantCulture));
                        Tmean.Add(float.Parse(values[1], CultureInfo.InvariantCulture));
                        Tdew.Add(float.Parse(values[2], CultureInfo.InvariantCulture));
                        Tground.Add(float.Parse(values[3], CultureInfo.InvariantCulture));
                        RH.Add(float.Parse(values[4], CultureInfo.InvariantCulture));
                        AP.Add(float.Parse(values[5], CultureInfo.InvariantCulture));
                        u.Add(float.Parse(values[6], CultureInfo.InvariantCulture));
                        Rad.Add(ParseFloatWithNaN(values[7]));
                        Rad_MJm2d.Add(ParseFloatWithNaN(values[8]));
                        Sunshine.Add(ParseFloatWithNaN(values[9]));
                        PET.Add(float.Parse(values[10], CultureInfo.InvariantCulture));
                        Precipitation_mm_even.Add(float.Parse(values[11], CultureInfo.InvariantCulture));
                        Discharge.Add(ParseFloatWithNaN(values[12]));
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("Error parsing line: " + line);
                        Console.WriteLine(ex.Message);
                    }
                }
            }
        }

        private CsvMeteoData() { }

        static float ParseFloatWithNaN(string value)
        {
            return value.Equals("NaN", StringComparison.OrdinalIgnoreCase) ? float.NaN : float.Parse(value, CultureInfo.InvariantCulture);
        }

        public CsvMeteoData GetSubRange(int start, int length)
        {
            int endExclusive = start + length;
            var subRange = new CsvMeteoData();
            subRange.allDateTimes = allDateTimes.GetRange(start, endExclusive - start);
            subRange.Tmean = Tmean.GetRange(start, endExclusive - start);
            subRange.Tdew = Tdew.GetRange(start, endExclusive - start);
            subRange.Tground = Tground.GetRange(start, endExclusive - start);
            subRange.RH = RH.GetRange(start, endExclusive - start);
            subRange.AP = AP.GetRange(start, endExclusive - start);
            subRange.u = u.GetRange(start, endExclusive - start);
            subRange.Rad = Rad.GetRange(start, endExclusive - start);
            subRange.Rad_MJm2d = Rad_MJm2d.GetRange(start, endExclusive - start);
            subRange.Sunshine = Sunshine.GetRange(start, endExclusive - start);
            subRange.PET = PET.GetRange(start, endExclusive - start);
            subRange.Precipitation_mm_even = Precipitation_mm_even.GetRange(start, endExclusive - start);
            subRange.Discharge = Discharge.GetRange(start, endExclusive - start);
            return subRange;
        }
    }
}
