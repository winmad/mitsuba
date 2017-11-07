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

class Heightfield2Trimesh : public Utility {
public:
	int run(int argc, char **argv) {
		Properties props = Properties("bitmap");
		props.setString("filename", argv[1]);
		props.setString("wrapMode", "repeat");
		Texture *texture = static_cast<Texture *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Texture), props));
		texture->configure();

		props = Properties("heightfield");
		props.setTransform("toWorld", Transform::scale(Vector3(std::atof(argv[2]), std::atof(argv[3]), 1.0f)));
		props.setFloat("scale", 1.0f);
		Shape *heightfield = static_cast<Shape *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Shape), props));
		heightfield->addChild(texture);
		heightfield->configure();

		ref<TriMesh> trimesh = heightfield->createTriMesh();
		fs::path filename(argv[1]);
		
		//filename.replace_extension(".obj");
		//trimesh->writeOBJ(filename);
		
		filename.replace_extension(".trimesh");
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);
		trimesh->serialize(stream);
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(Heightfield2Trimesh, "Convert heightfield to trimesh")
MTS_NAMESPACE_END
