#include <mitsuba/render/util.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/volume.h>
#include <iomanip>
#include "fluence2_proc.h"


MTS_NAMESPACE_BEGIN


class FluenceTracer2: public Utility
{
public:
    FluenceTracer2()
    {
    }


    static void saveEXR(const Spectrum *img, Vector2i size, const char *f)
    {
        ref<Bitmap> bmp = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, size);
        Vector4* data = reinterpret_cast<Vector4*>(bmp->getFloatData());
        int index = 0;
        for ( int i = 0; i < size.x*size.y; i++ )
        {
            const Spectrum &s = img[index];
            s.toLinearRGB(data[index][0], data[index][1], data[index][2]);
            data[index++][3] = 1.0f;
            //data[index++] = Vector4(s[0], s[1], s[2], 1.0f);
        }

        ref<FileStream> stream = new FileStream(f, FileStream::ETruncWrite);
        bmp->write(Bitmap::EOpenEXR, stream);
    }


    static std::string matrixToString(const Matrix4x4 &mat)
    {
        std::ostringstream oss;
        for ( int i = 0; i < 4; ++i )
            for ( int j = 0; j < 4; ++j )
            {
                oss.setf(std::ios::fixed);
                oss.precision(6);
                oss << mat(i, j) << ' ';
            }
        return oss.str();
    }


    static Float getScale(const Matrix4x4 &mat)
    {
        return mat.col(0).length()*mat.col(1).length()*mat.col(2).length();
    }


    bool checkScene(const Scene *scene, AABB &aabb)
    {
        ref_vector<Shape> shapes = scene->getShapes();
        size_t nshapes = shapes.size();
        size_t nblock = 0, nemitter = 0;
        for ( size_t i = 0; i < nshapes; ++i )
            if ( shapes[i]->getInteriorMedium() )
            {
                ++nblock;
                aabb = shapes[i]->getAABB();
            }
            else if ( shapes[i]->getEmitter() )
                ++nemitter;
        if ( nblock != 1 )
        {
            std::cout << "Cannot determine the medium block." << std::endl;
            return false;
        }
        else
            std::cout << "Medium " << aabb.toString() << std::endl;

        ref_vector<ConfigurableObject> refobjs = scene->getReferencedObjects();
        if ( refobjs.size() != 1 || !refobjs[0]->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource)) )
        {
            std::cout << "Failed to detect the density volume. " <<
                "Need to declare it directly under <Scene> and refer to it in <Medium>." << std::endl;
            return false;
        }

        return true;
    }


    void generateVPLTasks(const VolumeDataSource *densityVol, const AABB &aabb, const Vector3u &reso,
        std::vector<Vector3u> &buf_task, std::vector<Float> &buf_weight)
    {
        Vector aabbExtents = aabb.getExtents();
        ref<Sampler> sampler0 = static_cast<Sampler *>(PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("independent")));

        const int nsample = 5000;
        buf_task.clear();
        buf_weight.clear();
        for ( uint32_t z = 0; z < reso.z; ++z )
            for ( uint32_t y = 0; y < reso.y; ++y )
                for ( uint32_t x = 0; x < reso.x; ++x )
                {
                    int tot = 0;
                    for ( int i = 0; i < nsample; ++i )
                    {
                        Point p;
                        p.x = aabb.min.x + (static_cast<Float>(x) + sampler0->next1D())*aabbExtents.x
                            /static_cast<Float>(reso.x);
                        p.y = aabb.min.y + (static_cast<Float>(y) + sampler0->next1D())*aabbExtents.y
                            /static_cast<Float>(reso.y);
                        p.z = aabb.min.z + (static_cast<Float>(z) + sampler0->next1D())*aabbExtents.z
                            /static_cast<Float>(reso.z);

                        if ( densityVol->lookupFloat(p) > 0.1f )
                            ++tot;
                    }

                    // XXX: make 5% changeable
                    if ( static_cast<Float>(tot) > 0.05f*static_cast<Float>(nsample) )
                        buf_task.push_back(Vector3u(x, y, z));
                    buf_weight.push_back(
                        tot ? static_cast<Float>(nsample)/static_cast<Float>(tot) : 0.0f
                    );
                }
    }


    int run(int argc, char **argv)
    {
        size_t nphoton;
        std::string outputDir = "out";
        size_t granularity = 10000;
        
        if ( argc < 6 || (nphoton = atoi(argv[2])) <= 0 )
        {
            std::cout << "Usage: mtsutil fluence2 [scene file] [number of photons]\n"
                << "[resoSurf] [resoVol] [resoVPL] [output directory] [granularity]" << std::endl;
            return 1;
        }

        int iResoSurf = atoi(argv[3]);
        int iResoVol = atoi(argv[4]);
        int iResoVPL = atoi(argv[5]);
        if ( iResoSurf <= 0 || iResoVol <= 0 || iResoVPL <= 0 )
        {
            std::cout << "Invalid resolution" << std::endl;
            return 2;
        }

        if ( argc > 6 ) outputDir = argv[6];
        std::cout << "Output Directory: [" << outputDir << ']' << std::endl;

        if ( argc > 7 ) granularity = static_cast<size_t>(atoi(argv[7]));
        std::cout << "Granularity: " << granularity << std::endl;

        fs::path scene_file;
        const int maxDepth = -1;
        Vector3u resoSurf, resoVol, resoVPL;
        AABB mediumAABB;
        const VolumeDataSource *densityVol;
        std::vector<Vector3u> buf_task;
        std::vector<Float> buf_weight;
        uint32_t ntask;
        ref<Scene> scene0;

        {
            std::string fn0 = argv[1];
            ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
            scene_file = fileResolver->resolve(fn0);
            Utility::ParameterMap paramMap;
            paramMap["mat"] = matrixToString(
                Transform::translate(Vector(0.0f, 0.0f, 10000.0f)).getMatrix()
            );
            scene0 = loadScene(scene_file.string(), paramMap);
            scene0->incRef();
            scene0->initialize();

            if ( !checkScene(scene0.get(), mediumAABB) ) return 2;
            
            // synthetic
            //resoSurf = CrossImage::computeResolution(mediumAABB, 10);
            //resoVol = CrossImage::computeResolution(mediumAABB, 16);
            //resoVPL = CrossImage::computeResolution(mediumAABB, 16);

            // felt
            //resoSurf = CrossImage::computeResolution(mediumAABB, 25);
            //resoVol = CrossImage::computeResolution(mediumAABB, 50);
            //resoVPL = CrossImage::computeResolution(mediumAABB, 50);

            // gabardine
            //resoSurf = CrossImage::computeResolution(mediumAABB, 15);
            //resoVol = CrossImage::computeResolution(mediumAABB, 25);
            //resoVPL = CrossImage::computeResolution(mediumAABB, 25);

            resoSurf = CrossImage::computeResolution(mediumAABB, iResoSurf);
            resoVol = CrossImage::computeResolution(mediumAABB, iResoVol);
            resoVPL = CrossImage::computeResolution(mediumAABB, iResoVPL);

            std::cout << "resoSurf: " << resoSurf.toString() << '\n'
                << "resoVol: " << resoVol.toString() << '\n'
                << "resoVPL: " << resoVPL.toString() << std::endl;

            densityVol = static_cast<const VolumeDataSource *>(scene0->getReferencedObjects()[0].get());
            generateVPLTasks(densityVol, mediumAABB, resoVPL, buf_task, buf_weight);
            ntask = static_cast<uint32_t>(buf_task.size());
        }

        // Save the task list to a file
        {
            FILE *fout = fopen((outputDir + "/info.txt").c_str(), "wt");
            fprintf(fout, "%u %u %u\n", resoSurf.x, resoSurf.y, resoSurf.z);
            fprintf(fout, "%u %u %u\n", resoVol.x, resoVol.y, resoVol.z);
            fprintf(fout, "%u %u %u\n", resoVPL.x, resoVPL.y, resoVPL.z);
            fprintf(fout, "%u\n", ntask);
            for ( uint32_t i = 0; i < ntask; ++i )
                fprintf(fout, "%u %u %u\n", buf_task[i].x, buf_task[i].y, buf_task[i].z);
            fclose(fout);

            ref<FileStream> stream = new FileStream(outputDir + "/weight.bin", FileStream::ETruncWrite);
            stream->writeSize(buf_weight.size());
            stream->writeFloatArray(&buf_weight[0], buf_weight.size());
        }

        CrossImage crossImg(mediumAABB, resoSurf);
        uint32_t maxIdx = crossImg.getMaxIndex();

        //Spectrum *transferMat = new Spectrum[maxIdx*maxIdx];
        //memset(transferMat, 0, maxIdx*maxIdx*sizeof(Spectrum));
        std::vector<Spectrum> transSS, transSV;
        std::vector<Spectrum> transVS, transVV;
        transSS.resize(maxIdx*maxIdx);
        transSV.resize(maxIdx*resoVol.x*resoVol.y*resoVol.z);
        transVS.resize(ntask*maxIdx);
        transVV.resize(ntask*resoVol.x*resoVol.y*resoVol.z);

        std::string fnameSS = outputDir + "/transSS.exr";
        std::string fnameSV = outputDir + "/transSV.exr";
        std::string fnameVS = outputDir + "/transVS.exr";
        std::string fnameVV = outputDir + "/transVV.exr";

        for ( uint32_t idx = 0; idx < ntask + maxIdx; ++idx )
        {
            
            Float emitterArea = 0.0f;
            AABB emitterAABB;
            ref<Scene> scene;

            if ( idx < ntask )
            {
                scene = scene0;
                Vector aabbExtents = mediumAABB.getExtents();
                emitterAABB.min.x = mediumAABB.min.x + static_cast<Float>(buf_task[idx].x)*aabbExtents.x
                    /static_cast<Float>(resoVPL.x);
                emitterAABB.min.y = mediumAABB.min.y + static_cast<Float>(buf_task[idx].y)*aabbExtents.y
                    /static_cast<Float>(resoVPL.y);
                emitterAABB.min.z = mediumAABB.min.z + static_cast<Float>(buf_task[idx].z)*aabbExtents.z
                    /static_cast<Float>(resoVPL.z);
                emitterAABB.max.x = emitterAABB.min.x + aabbExtents.x/static_cast<Float>(resoVPL.x);
                emitterAABB.max.y = emitterAABB.min.y + aabbExtents.y/static_cast<Float>(resoVPL.y);
                emitterAABB.max.z = emitterAABB.min.z + aabbExtents.z/static_cast<Float>(resoVPL.z);
            }
            else
            {
                const Transform trans = crossImg.indexToTransform(idx - ntask);
                Utility::ParameterMap paramMap;
                paramMap["mat"] = matrixToString(trans.getMatrix());
                emitterArea = 4.0f*getScale(trans.getMatrix());
                scene = loadScene(scene_file.string(), paramMap);
                scene->initialize();
            }

            if ( idx < ntask )
                std::cout << "Emitter " << emitterAABB.toString() << std::endl;
            else
                std::cout << std::setprecision(6) << "Emitter surface area: " << emitterArea << std::endl;

            ref<Scheduler> sched = Scheduler::getInstance();
            int sceneResID = sched->registerResource(scene);
            int sensorResID = sched->registerResource(scene->getSensor());

            ref<Sampler> sampler = static_cast<Sampler *>(PluginManager::getInstance()->
                createObject(MTS_CLASS(Sampler), Properties("independent")));
            //int samplerResID = sched->registerResource(sampler);
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

            std::ostringstream oss;
            if ( idx < ntask )
                oss << "VPL " << (idx + 1) << '/' << ntask;
            else
                oss << "Boundary " << (idx - ntask + 1) << '/' << maxIdx;

            ref<Timer> timer = new Timer();
            timer->reset();
            ref<Flunence2ParticleProcess> proc = idx < ntask
                ? new Flunence2ParticleProcess(nphoton, granularity, maxDepth, densityVol, mediumAABB,
                    resoSurf, resoVol, emitterAABB, oss.str().c_str())
                : new Flunence2ParticleProcess(nphoton, granularity, maxDepth, densityVol, mediumAABB,
                    resoSurf, resoVol, emitterArea, oss.str().c_str());
            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", samplerResID);
            scene->bindUsedResources(proc);

            sched->schedule(proc);
            sched->wait(proc);

            std::cout << "\nFinished in " << static_cast<Float>(timer->getMilliseconds()/1000.0)
                << " secs." << std::endl;

            const BlockInfo& out = proc->getOutput();
#if 0
            char fname0[255] = {0}, fname1[255] = {0};
            if ( idx < ntask )
            {
                sprintf(fname0, "%s/volume_%04u.exr", outputDir.c_str(), idx);
                sprintf(fname1, "%s/volume_%04u.vol", outputDir.c_str(), idx);
            }
            else
            {
                sprintf(fname0, "%s/boundary_%04u.exr", outputDir.c_str(), idx - ntask);
                sprintf(fname1, "%s/boundary_%04u.vol", outputDir.c_str(), idx - ntask);
            }
            out.saveCrossImage(fname0);
            out.saveVOL(fname1);
#endif
            const Spectrum *rawSurf = out.getRawSurfaceData();
            const Spectrum *rawVol = out.getRawVolumeData();
            if ( idx < ntask )
            {
                for ( uint32_t i = 0; i < maxIdx; ++i )
                    transVS[i*ntask + idx] = rawSurf[i];
                for ( uint32_t i = 0; i < resoVol.x*resoVol.y*resoVol.z; ++i )
                    transVV[i*ntask + idx] = rawVol[i];

                if ( idx % 20 == 0 )
                {
                    saveEXR(&transVS[0], Vector2i(ntask, maxIdx), fnameVS.c_str());
                    saveEXR(&transVV[0], Vector2i(ntask, resoVol.x*resoVol.y*resoVol.z), fnameVV.c_str());
                }
            }
            else
            {
                for ( uint32_t i = 0; i < maxIdx; ++i )
                    transSS[i*maxIdx + idx - ntask] = rawSurf[i];
                for ( uint32_t i = 0; i < resoVol.x*resoVol.y*resoVol.z; ++i )
                    transSV[i*maxIdx + idx - ntask] = rawVol[i];

                if ( (idx - ntask) % 20 == 0 )
                {
                    saveEXR(&transSS[0], Vector2i(maxIdx, maxIdx), fnameSS.c_str());
                    saveEXR(&transSV[0], Vector2i(maxIdx, resoVol.x*resoVol.y*resoVol.z), fnameSV.c_str());
                }
            }
        }

        //saveEXR(transferMat, Vector2i(maxIdx), "out\\transfer.exr");
        //delete[] transferMat;
        saveEXR(&transSS[0], Vector2i(maxIdx, maxIdx), fnameSS.c_str());
        saveEXR(&transSV[0], Vector2i(maxIdx, resoVol.x*resoVol.y*resoVol.z), fnameSV.c_str());
        saveEXR(&transVS[0], Vector2i(ntask, maxIdx), fnameVS.c_str());
        saveEXR(&transVV[0], Vector2i(ntask, resoVol.x*resoVol.y*resoVol.z), fnameVV.c_str());

        scene0->decRef();

        return 0;
    }

    MTS_DECLARE_UTILITY()
};


MTS_EXPORT_UTILITY(FluenceTracer2, "Fluence particle tracer")
MTS_NAMESPACE_END
