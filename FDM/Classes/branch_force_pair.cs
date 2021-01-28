using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace FDM.Classes
{
    class branch_force_pair
    {
        // Member variables
        static public int count = 0;
        public int id { get; set; }
        public Line line { get; set; }
        public double ratio { get; set; }

        // Constructor
        public branch_force_pair(Line _line, double _ratio)
        {
            id = count;
            line = _line;
            ratio = _ratio;
            count += 1;
        }
    }
}
