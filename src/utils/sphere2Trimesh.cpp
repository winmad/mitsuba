#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <boost/filesystem/path.hpp>

MTS_NAMESPACE_BEGIN

class Sphere2Trimesh : public Utility {
public:
	int run(int argc, char **argv) {
		Properties props;

		props = Properties("sphere");
		props.setFloat("radius", 1.0);
		Shape *sphere = static_cast<Shape *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Shape), props));
		sphere->configure();

		ref<TriMesh> trimesh = sphere->createTriMesh();
		fs::path filename(argv[1]);

		filename.replace_extension(".obj");
		trimesh->writeOBJ(filename);

		//filename.replace_extension(".trimesh");
		//ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		//trimesh->serialize(stream);
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(Sphere2Trimesh, "Convert sphere to trimesh")
MTS_NAMESPACE_END
