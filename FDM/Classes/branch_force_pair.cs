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
        public int id { get; set; }
        public Line line { get; set; }
        public double ratio { get; set; }

        // Constructor

        public branch_force_pair(Line _line, double _ratio, int _id)
        {
            id = _id;
            line = _line;
            ratio = _ratio;
        }
    }
}
