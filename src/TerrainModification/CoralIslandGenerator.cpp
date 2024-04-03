#include "CoralIslandGenerator.h"

#include "DataStructure/Matrix3.h"
#include "Graphics/DisplayGraphics.h"

CoralIslandGenerator::CoralIslandGenerator()
{

}

GridF CoralIslandGenerator::generate(GridF heights, float subsidence, float waterLevel, float minCoralHeight,
                                              float maxCoralHeight, float verticalScale, float horizontalScale, float alpha)
{
    GridF perlinMap = GridF::perlin(heights.getDimensions(), Vector3(10.f, 10.f, 10.f)).normalize();

    float vh = 1.0 / horizontalScale;
    float vv = 1.0 / verticalScale;
//    float alpha = 0.1;

//    GridF heights = heightmap.getHeights();
    float downscale = heights.max();
    heights /= downscale;

    GridF initialCoral(heights.getDimensions());
    GridF lowerBand(heights.getDimensions());

    initialCoral = heights.binarizeBetween(minCoralHeight, maxCoralHeight, true, true);
    GridF dilatedInitialSeed = initialCoral.dilate();

    for (size_t i = 0; i < lowerBand.size(); i++)
        lowerBand[i] = (initialCoral[i] != dilatedInitialSeed[i] && heights[i] <= minCoralHeight ? 1.f : lowerBand[i]);
    GridF distanceFromLowerCorals = (1.f - lowerBand).toDistanceMap(true, false); //.normalized();
    distanceFromLowerCorals = distanceFromLowerCorals.max() - distanceFromLowerCorals;
    for (size_t i = 0; i < distanceFromLowerCorals.size(); i++) {
        if (heights[i] > minCoralHeight)
            distanceFromLowerCorals[i] *= -1.f;
    }

    float distanceWhenLowerCoralTouchesWater = (vh / vv) * (maxCoralHeight - minCoralHeight); // / verticalScale;

    GridI heightAboveWater = heights.binarize(waterLevel, true, true);
    GridI distanceIsSmall = distanceFromLowerCorals.binarize(distanceWhenLowerCoralTouchesWater, false, true);
    GridI heightLowerThanMaxCoral = heights.binarize(maxCoralHeight, false, true);
    GridI heightAboveLowBand = heights.binarize(minCoralHeight, true, true);
    GridI heightAboveHighBand = heights.binarize(maxCoralHeight, true, true);


    GridF insideCorals(heights.getDimensions());
    for (size_t i = 0; i < insideCorals.size(); i++)
        if (heightAboveLowBand[i]) //if (heightAboveWater[i] || (distanceIsSmall[i] && !heightAboveLowBand[i]))
            insideCorals[i] = 1.f;
    GridF aFactorIHaveToCheck = maxCoralHeight - (1.f - (distanceFromLowerCorals.abs() / (distanceFromLowerCorals.abs() * insideCorals).max())) * (maxCoralHeight - minCoralHeight) * alpha;

    for (size_t i = 0; i < insideCorals.size(); i++)
        insideCorals[i] *= aFactorIHaveToCheck[i];

    GridF outsideCorals(heights.getDimensions());
    for (size_t i = 0; i < outsideCorals.size(); i++)
        if (!heightAboveLowBand[i]) //if (distanceFromLowerCorals[i] >= distanceWhenLowerCoralTouchesWater)
            outsideCorals[i] = 1.f;
    GridF aFactorIHaveToCheck2 = (maxCoralHeight - (1.f - (distanceFromLowerCorals.abs() / (distanceFromLowerCorals.abs() * insideCorals).max()) * alpha));
//    GridF aFactorIHaveToCheck2 = waterLevel + (distanceWhenLowerCoralTouchesWater - distanceFromLowerCorals) / (distanceWhenLowerCoralTouchesWater * 5.f);
    for (size_t i = 0; i < insideCorals.size(); i++)
        outsideCorals[i] *= aFactorIHaveToCheck2[i];

    GridF finalMap = GridF::max(GridF::max(insideCorals, outsideCorals)/* * clamp(std::pow(1.f - subsidence, .5f) + .8f, 0.f, 1.f)*/, heights * subsidence);
    return finalMap * downscale;
}

std::vector<EnvObject*> CoralIslandGenerator::envObjsFromFeatureMap(const GridV3& img, const Vector3& terrainDimensions)
{
    // Input image might not be the same size than the terrain, need to resize all the curves on the XY components
    Vector3 ratio = terrainDimensions / img.getDimensions();
    ratio.z = 1;

    // Map the image color to a type of object
    std::map<std::tuple<int, int, int>, std::string> colorToFeature = {
        {{255,   0,   0}, "abyss"},
        {{  0,   0, 255}, "reef"},
        {{  0, 255, 255}, "lagoon"},
        {{  0, 255,   0}, "coast"},
        {{255, 255,   0}, "island"}
    };

    std::map<std::string, GridI> featureAreas;

    displayProcessTime("Extracting shapes... ", [&]() {
        // Create binary masks for each of the objects
        for (auto& [_, name] : colorToFeature)
            featureAreas[name] = GridI(img.getDimensions());

        for (size_t i = 0; i < img.size(); i++) {
            const auto& pix = img[i];
            if (colorToFeature.count({pix.x, pix.y, pix.z}) == 0) continue;
            featureAreas[colorToFeature[{pix.x, pix.y, pix.z}]][i] = 1;
        }
    });

    displayProcessTime("Filling holes in shapes... ", [&]() {
        // If some binary masks have holes (e.g. the island is inside the lagoon), remove them.
        // This is done using CCL algorithm, maybe not the fastest and smartest way
        for (auto& [name, area] : featureAreas) {
            if (name == "reef") continue;
            area = area.fillHoles(true);
        }
    });


    std::vector<EnvObject*> objects;
/*
    std::cout << "Reefs: " << showTime(timeIt([&]() {
        // Extract the lagoon contours to instantiate the lagoons and the reefs
        auto skeletons = featureAreas["reef"].skeletonizeToBSplines();
        for (auto& curve : skeletons) {
            curve.scale(ratio);
            BSpline simplifiedCurve = curve;
            simplifiedCurve = simplifiedCurve.getPath(50); // Reduce the complexity of the curve to avoid having too much computations after
            if (simplifiedCurve.length() < 5.f) continue; // Remove too small elements

            EnvCurve* reef = dynamic_cast<EnvCurve*>(EnvObject::instantiate("reef"));
            reef->curve = simplifiedCurve;
            objects.push_back(reef);
        }
    })) << std::endl;
*/

    displayProcessTime("Finding lagoons and reefs... ", [&]() {
        // Extract the lagoon contours to instantiate the lagoons and the reefs
        int nbLagoonsKept = 0; // Keep track of the number of lagoons reaching threshold area
        auto lagoonContours = featureAreas["lagoon"].findContoursAsCurves();
        for (auto& curve : lagoonContours) {
            curve.scale(ratio);
            ShapeCurve simplifiedCurve = curve;
            simplifiedCurve = simplifiedCurve.getPath(10); // Reduce the complexity of the curve to avoid having too much computations after
            if (simplifiedCurve.computeArea() < 15.f) continue; // Remove too small elements

            EnvCurve* reef = dynamic_cast<EnvCurve*>(EnvObject::instantiate("reef"));
            reef->curve = simplifiedCurve;
            objects.push_back(reef);

            EnvArea* lagoon = dynamic_cast<EnvArea*>(EnvObject::instantiate("lagoon"));
            lagoon->area = simplifiedCurve;
            objects.push_back(lagoon);

            nbLagoonsKept++;
        }

        if (nbLagoonsKept == 0) {  // In the case where we have reefs but no lagoon
            auto skeletons = featureAreas["reef"].skeletonizeToBSplines();
            std::cout << "Skeletons : " << skeletons.size() << std::endl;
            for (auto& curve : skeletons) {
                curve.scale(ratio);
                BSpline simplifiedCurve = curve;
                simplifiedCurve = simplifiedCurve.getPath(50); // Reduce the complexity of the curve to avoid having too much computations after
                if (simplifiedCurve.length() < 5.f) continue; // Remove too small elements

//                EnvCurve* reef = dynamic_cast<EnvCurve*>(EnvObject::instantiate("reef"));
//                reef->curve = simplifiedCurve;
                EnvArea* reef = dynamic_cast<EnvArea*>(EnvObject::instantiate("reef"));
                reef->area = simplifiedCurve;
                objects.push_back(reef);
            }
        }
    });

    displayProcessTime("Finding coasts... ", [&]() {
        // Extract the coast contours
        auto coastContours = featureAreas["coast"].findContoursAsCurves();
        for (auto& curve : coastContours) {
            curve.scale(ratio);
            ShapeCurve simplifiedCurve = curve;
            simplifiedCurve = simplifiedCurve.getPath(50); // Reduce the complexity of the curve to avoid having too much computations after
            if (simplifiedCurve.computeArea() < 15.f) continue; // Remove too small elements

            EnvArea* coast = dynamic_cast<EnvArea*>(EnvObject::instantiate("coast"));
            coast->area = simplifiedCurve;
            objects.push_back(coast);
        }
    });


    displayProcessTime("Finding islands... ", [&]() {
        // Extract the island contours
        auto islandContours = featureAreas["island"].findContoursAsCurves();
        for (auto& curve : islandContours) {
            curve.scale(ratio);
            ShapeCurve simplifiedCurve = curve;
            simplifiedCurve = simplifiedCurve.getPath(20); // Reduce the complexity of the curve to avoid having too much computations after
            if (simplifiedCurve.computeArea() < 15.f) continue; // Remove too small elements
            EnvArea* island = dynamic_cast<EnvArea*>(EnvObject::instantiate("island"));
            island->area = simplifiedCurve;
            objects.push_back(island);
        }
    });

    std::cout << objects.size() << " objects found" << std::endl;
    return objects;
}
