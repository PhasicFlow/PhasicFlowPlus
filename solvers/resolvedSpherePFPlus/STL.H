#include <sys/stat.h>

triSurface generateSphere(const vector& center, const scalar& diameter, int resolution)
{
    scalar radius = diameter / 2.0;
    pointField points;
		triFaceList faces;

    // Generate vertices
    for (int i = 0; i <= resolution; ++i) {
        scalar theta = M_PI * i / resolution;
        for (int j = 0; j <= resolution; ++j) {
            scalar phi = 2.0 * M_PI * j / resolution;
            scalar x = center.x() + radius * Foam::sin(theta) * Foam::cos(phi);
            scalar y = center.y() + radius * Foam::sin(theta) * Foam::sin(phi);
            scalar z = center.z() + radius * Foam::cos(theta);
            points.append(point(x, y, z));
        }
    }

    // Generate faces
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            int p1 = i * (resolution + 1) + j;
            int p2 = p1 + 1;
            int p3 = p1 + (resolution + 1);
            int p4 = p3 + 1;

						// To generate starting triangles
						if (i == 0)
						{
                faces.append(triFace(labelList({p1, p3, p4})));
            }
						// To generate ending triangles
						else if (i == resolution - 1)
						{
								faces.append(triFace(labelList({p1, p4, p2})));
						}
						// To generate the rest
						else
						{
                faces.append(triFace(labelList({p1, p3, p4})));
								faces.append(triFace(labelList({p1, p4, p2})));
            }
        }
    }

    // Create the triSurface
		triSurface sphere(faces, points);
    sphere.cleanup(false);

		return sphere;
}

void writeSTLs(PtrList<triSurface>& STLs, word timeName)
{
		// Checking if the STL folder exist
		struct stat info;
		const std::string path = "STL/";

		if (stat(path.c_str(), &info) != 0)
		{
				mkdir(path.c_str(), 0755);
		}

		OFstream outputFile("STL/particles_" + timeName + ".stl");

		forAll (STLs, i)
		{
				triSurface STLi(STLs[i]);

	      outputFile  << "solid " << "particle" + std::to_string(i) << endl;

	      forAll(STLi, facei)
	      {
	          const vector& n = STLi[facei].unitNormal(STLi.points());

	          outputFile  << "  facet normal "
	              << n.x() << ' ' << n.y() << ' ' << n.z() << nl
	              << "    outer loop" << endl;

	          point pa = STLi[facei].points(STLi.points())[0];
						point pb = STLi[facei].points(STLi.points())[1];
						point pc = STLi[facei].points(STLi.points())[2];

	          outputFile  << "       vertex "
	              << pa.x() << ' ' << pa.y() << ' ' << pa.z() << nl
	              << "       vertex "
	              << pb.x() << ' ' << pb.y() << ' ' << pb.z() << nl
	              << "       vertex "
	              << pc.x() << ' ' << pc.y() << ' ' << pc.z() << nl
	              << "    endloop" << nl
	              << "  endfacet" << endl;
	      }

	      outputFile  << "endsolid " << "particle" + std::to_string(i) << endl;
		}
}
