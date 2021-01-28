using System;
using System.Collections.Generic;
using FDM.Classes;
using Grasshopper.Kernel;
using Rhino.Geometry;
using la = MathNet.Numerics.LinearAlgebra;


// In order to load the result of this wizard, you will also need to
// add the output bin/ folder of this project to the list of loaded
// folder in Grasshopper.
// You can use the _GrasshopperDeveloperSettings Rhino command for that.

namespace FDM
{
    public class FDMComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public FDMComponent()
          : base("Force Density Method", "FDM",
              "Force Density Method",
              "Form Finding", "FDM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddLineParameter("Lines", "l", "Lines representing individual bars in grid", GH_ParamAccess.list);
            pManager.AddNumberParameter("Force densities", "q", "The force densities of the bars", GH_ParamAccess.list);
            pManager.AddPointParameter("Support points", "sp", "The strucutres suport points", GH_ParamAccess.list);
            pManager.AddVectorParameter("Force Vector", "F", "The force applied to the nodes", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Lines", "L", "The output lines from FDM", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Retrieving input params
            List<Line> lines = new List<Line>();
            List<double> qs = new List<double>();
            List<Point3d> sPts = new List<Point3d>();
            Vector3d fVec = new Vector3d();

            if (!DA.GetDataList(0, lines)) return;
            if (!DA.GetDataList(1, qs)) return;
            if (!DA.GetDataList(2, sPts)) return;
            if (!DA.GetData(3, ref fVec)) return;

            List<branch_force_pair> bf_pairs = new List<branch_force_pair>();
            for (int i = 0; i < lines.Count; i++)
            {
                bf_pairs.Add(new branch_force_pair(lines[i], qs[i], i));
            }

            List<Point3d> sortedPoints = SortPoints(bf_pairs, sPts);
            la.Matrix<double> cMat = BuildBranchNodeMatrix(bf_pairs, sortedPoints);
            la.Vector<double> qVec = ForceDensityVector(bf_pairs);

            List<la.Matrix<double>> matrices = dMatrices(cMat, qVec, sPts);
            la.Matrix<double> dN = matrices[0];
            la.Matrix<double> dF = matrices[1];

            List<la.Vector<double>> coordVec = CoordinateVectors(sortedPoints);


            la.Vector<double> xf = coordVec[0].SubVector(coordVec[0].Count - sPts.Count, sPts.Count);
            la.Vector<double> yf = coordVec[1].SubVector(coordVec[1].Count - sPts.Count, sPts.Count);
            la.Vector<double> zf = coordVec[2].SubVector(coordVec[2].Count - sPts.Count, sPts.Count);

            la.Vector<double> px = la.Double.Vector.Build.Dense(sortedPoints.Count - sPts.Count, fVec.X);
            la.Vector<double> py = la.Double.Vector.Build.Dense(sortedPoints.Count - sPts.Count, fVec.Y);
            la.Vector<double> pz = la.Double.Vector.Build.Dense(sortedPoints.Count - sPts.Count, fVec.Z);

            la.Vector<double> xN = dN.Solve(px - dF * xf);
            la.Vector<double> yN = dN.Solve(py - dF * yf);
            la.Vector<double> zN = dN.Solve(pz - dF * zf);

            List<Line> new_lines = newLines(cMat, xN, yN, zN, sPts);

            // Set the data
            DA.SetDataList(0, new_lines);

        }

        private List<Point3d> SortPoints(List<branch_force_pair> Branches, List<Point3d> supPoints)
        {
            List<Point3d> sortedPts = new List<Point3d>();
            foreach (Point3d p in supPoints) { sortedPts.Add(p); }

            foreach (branch_force_pair b in Branches)
            {
                Point3d spt = b.line.PointAt(0);
                Point3d ept = b.line.PointAt(1);
                bool addit1 = true;
                bool addit2 = true;
                foreach (Point3d p in sortedPts)
                {
                    if (p.DistanceTo(spt) < 0.001) addit1 = false;
                    if (p.DistanceTo(ept) < 0.001) addit2 = false;
                }
                if (addit1) sortedPts.Add(spt);
                if (addit2) sortedPts.Add(ept);
            }
            sortedPts.RemoveRange(0, supPoints.Count);
            sortedPts.AddRange(supPoints);
            return sortedPts;
        }

        private la.Matrix<double> BuildBranchNodeMatrix(List<branch_force_pair> branches, List<Point3d> _sortedPoints)
        {
            la.Matrix<double> c = la.Double.Matrix.Build.Dense(branches.Count, _sortedPoints.Count, 0);
            foreach (branch_force_pair b in branches)
            {
                int row = b.id;
                Point3d sPt = b.line.PointAt(0);
                Point3d ePt = b.line.PointAt(1);
                int col1 = 0;
                int col2 = 0;
                int i = 0;
                foreach (Point3d p in _sortedPoints)
                {
                    if (sPt.DistanceTo(p) < 0.001) col1 = i;
                    if (ePt.DistanceTo(p) < 0.001) col2 = i;
                    i += 1;
                }
                c[row, col1] = 1;
                c[row, col2] = -1;
            }
            return c;
        }

        private la.Vector<double> ForceDensityVector(List<branch_force_pair> branches)
        {
            la.Vector<double> q = la.Double.Vector.Build.Dense(branches.Count);
            for (int i = 0; i < branches.Count; i++)
            {
                q[i] = branches[i].ratio;
            }

            return q;
        }
        private List<la.Vector<double>> CoordinateVectors(List<Point3d> sortedPts)
        {
            la.VectorBuilder<double> v = la.Double.Vector.Build;
            la.Vector<double> x = v.Dense(sortedPts.Count);
            la.Vector<double> y = v.Dense(sortedPts.Count);
            la.Vector<double> z = v.Dense(sortedPts.Count);
            for (int i = 0; i < sortedPts.Count; i++)
            {
                x[i] = sortedPts[i].X;
                y[i] = sortedPts[i].Y;
                z[i] = sortedPts[i].Z;
            }
            List<la.Vector<double>> cordVecs = new List<la.Vector<double>>();
            cordVecs.Add(x);
            cordVecs.Add(y);
            cordVecs.Add(z);
            return cordVecs;
        }

        private List<la.Matrix<double>> dMatrices(la.Matrix<double> cMat, la.Vector<double> qVec, List<Point3d> supPts) {
            la.MatrixBuilder<double> m = la.Double.Matrix.Build;
            la.Matrix<double> qMat = m.DiagonalOfDiagonalVector(qVec);

            //Subtract the cN and cF of the C-matrix
            //la.Matrix<double> cN = la.Double.DenseMatrix.SubMatrix()
            la.Matrix<double> cN = cMat.SubMatrix(0, cMat.RowCount, 0, cMat.ColumnCount - supPts.Count);
            la.Matrix<double> cF = cMat.SubMatrix(0, cMat.RowCount, cMat.ColumnCount - supPts.Count, supPts.Count);

            la.Matrix<double> cNt = cN.Transpose();
            la.Matrix<double> dN = cNt * qMat * cN;
            
            la.Matrix<double> dF = cNt * qMat * cF;


            List<la.Matrix<double>> return_list = new List<la.Matrix<Double>>();
            return_list.Add(dN);
            return_list.Add(dF);
            // Return dN and dF as list -> return_list[0] is dN and return_list[1] is dF
            return return_list;
        }

        private List<Line> newLines(la.Matrix<double> cMat, la.Vector<double> xN, la.Vector<double> yN, la.Vector<double> zN, List<Point3d> supPts)
        {
            List<Point3d> newPts = new List<Point3d>();
            List<Line> newLines = new List<Line>();

            for (int i = 0; i < xN.Count; i++)
            {
                newPts.Add(new Point3d(xN[i], yN[i], zN[i]));
            }
            //Add the support points
            for (int i = 0; i < supPts.Count; i++)
            {
                newPts.Add(supPts[i]);
            }
            //Will this be correct?
            Point3d sPt = new Point3d();
            Point3d ePt = new Point3d();

            for (int i = 0; i < cMat.RowCount; i++)
            {
                for (int j = 0; j < cMat.ColumnCount; j++)
                {
                    if (cMat[i, j] == 1) { sPt = newPts[j];  }
                    if (cMat[i, j] == -1) { ePt = newPts[j]; }
                    
                }
            //Will sPt != ePt work in C#?
            if (sPt != ePt)
                {
                    newLines.Add(new Line(sPt, ePt));
                }
            }
            return newLines;
        }
        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("cdf3639a-f05a-4edd-8a73-e82832530970"); }
        }
    }
}
