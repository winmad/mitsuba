#include <mitsuba/render/util.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include "../../integrators/mft/fluencevol.h"
#include "../../volume/tetra2.h"
#include <Eigen/Dense>


MTS_NAMESPACE_BEGIN


class TestMFT: public Utility
{
public:
    TestMFT()
    {
    }


#if 0
    bool loadVOL(const std::string &fn, std::vector<Spectrum> &buf) const
    {
        ref<FileStream> stream = new FileStream(fn, FileStream::EReadOnly);

        // Header
        if ( stream->readChar() != 'V' ) return false;
        if ( stream->readChar() != 'O' ) return false;
        if ( stream->readChar() != 'L' ) return false;
        if ( stream->readUChar() != 0x3 ) return false;
        if ( stream->readInt() != 1 ) return false;

        // Volume size info.
        Vector3i reso(stream);
        if ( stream->readInt() != 3 ) return false;
        AABB aabb(stream);
        buf.resize(reso.x*reso.y*reso.z);

        // Volume data
        stream->readFloatArray(reinterpret_cast<Float *>(&buf[0]), buf.size()*3);

        return true;
    }
#endif


    void loadEXR(const std::string &fn, Eigen::MatrixXf &mat) const
    {
        ref<FileStream> stream = new FileStream(fn, FileStream::EReadOnly);
        ref<Bitmap> img = new Bitmap(Bitmap::EOpenEXR, stream);
        Assert( img->getPixelFormat() == Bitmap::ERGBA );
        const Vector4 *data = reinterpret_cast<const Vector4*>(img->getData());

        mat.resize(img->getHeight()*3, img->getWidth());
        int tot = 0;
        for ( int i = 0; i < img->getHeight(); ++i )
            for ( int j = 0; j < img->getWidth(); ++j )
            {
                mat(i, j) = data[tot][0];
                mat(img->getHeight() + i, j) = data[tot][1];
                mat(2*img->getHeight() + i, j) = data[tot][2];
                Assert( data[tot][3] == 1.0f );

                ++tot;
            }
    }


    int run(int argc, char **argv)
    {
        if ( argc < 4 )
        {
            std::cout << "Usage: mtsutil test_mft [scene file] [shell file] [info file]" << std::endl;
            return 1;
        }

#if 0
        Eigen::MatrixXf transSS, transSV, transVS, transVV;
        loadEXR("out4\\transSS.exr", transSS);
        loadEXR("out4\\transSV.exr", transSV);
        loadEXR("out4\\transVS.exr", transVS);
        loadEXR("out4\\transVV.exr", transVV);
#endif

        std::string fn = argv[1];
        ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
        fs::path scene_file = fileResolver->resolve(fn);
        ref<Scene> scene = loadScene(scene_file.string(), Utility::ParameterMap());
        scene->initialize();

        TetrahedronMesh mesh;
        if ( !mesh.load(argv[2]) )
        {
            std::cout << "Bad mesh" << std::endl;
            return 2;
        }
        mesh.configure();

        Properties fvolprop;
        fvolprop.setInteger("mode", 1);
        //fvolprop.setInteger("mode", 0);
        fvolprop.setString("filename", argv[3]);
        MFTFluenceVolume fvol(fvolprop);
        fvol.configure();

        Vector2u reso(fvol.getResolution());
        Vector3u blockReso(fvol.getVolumeResolution());

        uint32_t nblockX = reso.x;
        if ( argc > 4 )
        {
            nblockX = atoi(argv[4]);
            if ( nblockX > reso.x )
            {
                printf("Invalid nblockX value!\n");
                return 1;
            }
        }
        printf("nblockX = %u\n", nblockX);

        ref<Scheduler> sched = Scheduler::getInstance();
        int sceneResID = sched->registerResource(scene);
        int sensorResID = sched->registerResource(scene->getSensor());

        ref<Sampler> sampler = static_cast<Sampler *>(PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("independent")));
        std::vector<SerializableObject *> samplers(sched->getCoreCount());
        for ( size_t i = 0; i < sched->getCoreCount(); ++i )
        {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }
        int samplerResID = sched->registerMultiResource(samplers);
        for ( size_t i = 0; i < samplers.size(); ++i )
            samplers[i]->decRef();

        ref<MFTParticleProcess> proc = new MFTParticleProcess(
            50000000, 20000, 64, &mesh, reso, fvol.getVPLResolution()
        );
        proc->bindResource("scene", sceneResID);
        proc->bindResource("sensor", sensorResID);
        proc->bindResource("sampler", samplerResID);
        scene->bindUsedResources(proc);

        sched->schedule(proc);
        sched->wait(proc);
        std::cout << std::endl;

        const MFTParticleProcess::t_result &output = proc->getOutput();
        std::cout << "Stored " << proc->totalStoredParticles() << " particles." << std::endl;

        fvol.initialize(output);
        //fvol->propagate(0);

        ref<Bitmap> bmp = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(nblockX*blockReso.x, blockReso.z));

        AABB aabb = mesh.getAABB();
        std::cout << aabb.toString() << std::endl;
        Vector extents = aabb.getExtents();
        std::cout << aabb.toString() << std::endl;
        const int offsetX = (reso.x - nblockX)/2*blockReso.x;
        const int nsamples = 1000;
        const int niterations = 100;
        for ( int i = 0; i < static_cast<int>(blockReso.z); ++i )
        {
            std::cout << i << std::endl;
            for ( int j = 0; j < static_cast<int>(nblockX*blockReso.x); ++j )
            {
                Point p;
                p.x = aabb.min.x + extents.x*(static_cast<Float>(offsetX + j) + 0.5f)
                    /static_cast<Float>(reso.x*blockReso.x);
                p.y = aabb.min.y + extents.y*(static_cast<Float>(reso.y*blockReso.y/2) + 0.5f)
                    /static_cast<Float>(reso.y*blockReso.y);
                p.z = aabb.max.z - extents.z*(static_cast<Float>(i) + 0.5f)
                    /static_cast<Float>(blockReso.z);

                Spectrum s(0.0f);
                Point tex;
                if ( mesh.lookupPoint(p, tex) )
                    #pragma omp parallel for
                    for ( int k = 0; k < nsamples; ++k )
                    {
                        Spectrum val = fvol.lookupFluence(tex, sampler, &niterations);
                        #pragma omp critical
                        s += val;
                    }
                s *= static_cast<Float>(4.0*M_PI)/static_cast<Float>(nsamples);

                bmp->setPixel(Point2i(j, i), s);
            }
        }

        ref<FileStream> stream = new FileStream("out_mft.pfm", FileStream::ETruncWrite);
        bmp->write(Bitmap::EPFM, stream);

        std::cout << "done" << std::endl;

        return 0;
    }

    MTS_DECLARE_UTILITY()
};


MTS_EXPORT_UTILITY(TestMFT, "MFT Tester")
MTS_NAMESPACE_END
