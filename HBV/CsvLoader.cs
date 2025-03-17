using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HBV
{
    public class CsvLoader
    {
        public static List<IList<float>> LoadCsv(string filePath)
        {
            var result = new List<IList<float>>();

            foreach (var line in File.ReadLines(filePath))
            {
                //var values = line.Split(new char[] { ';', ',' }, 2,StringSplitOptions.RemoveEmptyEntries);
                var values = line.Split(',',StringSplitOptions.RemoveEmptyEntries);
                List<float> row = new List<float>();
                        for (int i = 0; i < values.Length; ++i)
                {
                    row.Add(Convert.ToSingle(values[i]));
                }         
                                 
                result.Add(row);
            }

            return result;
        }

        //static void Main(string[] args)
        //{
        //    string filePath = "data.csv"; // Change this to the actual file path
        //    List<List<double>> data = LoadCsv(filePath);

        //    // Example: Print loaded data
        //    foreach (var row in data)
        //    {
        //        Console.WriteLine(string.Join(", ", row));
        //    }
        //}
    }
}
