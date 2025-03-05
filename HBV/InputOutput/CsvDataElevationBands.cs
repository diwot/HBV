using System.Globalization;

namespace HBV.InputOutput
{
    public class CsvDataElevationBands
    {
        public List<float> ElevationBottom = new List<float>();
        public List<float> ElevationTop = new List<float>();
        public List<float> AreaM2 = new List<float>();
        public List<float> RelArea = new List<float>();
        public List<float> TotalAreaM2 = new List<float>();

        public CsvDataElevationBands(string path)
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
                        ElevationBottom.Add(float.Parse(values[0], CultureInfo.InvariantCulture));
                        ElevationTop.Add(float.Parse(values[1], CultureInfo.InvariantCulture));
                        AreaM2.Add(float.Parse(values[2], CultureInfo.InvariantCulture));
                        RelArea.Add(float.Parse(values[3], CultureInfo.InvariantCulture));
                        TotalAreaM2.Add(float.Parse(values[4], CultureInfo.InvariantCulture));
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("Error parsing line: " + line);
                        Console.WriteLine(ex.Message);
                    }
                }
            }
        }
    }
}
