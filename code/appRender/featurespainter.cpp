#include "featurespainter.h"
#include <QPainter>
#include <QPainterPath>
#include <iostream>

void remapTexCoordinates(const Vector2& p, const Box2& domain, const QSize& texSize, int& texi, int& texj)
{
    Vector2 q = p - domain.getMin();
    Vector2 d = domain.getMax() - domain.getMin();
    double u = q[0] / d[0];
    double v = q[1] / d[1];
    texi = int(u * (texSize.width() - 1));
    texj = int(v * (texSize.height() - 1));
}

double FeaturesPainter::randUniform(double a, double b)
{
    return a + (b-a)*uniformDist(rne);
}

QImage FeaturesPainter::baseTexture(const ScalarField2& hf, const ScalarField2& ice, const ScalarField2& ela) const
{
    double a, b;
    hf.getRange(a, b);
    Vector3 lightDir = Normalized(Vector3(-1.0, -1.0, 2.5));
    ColorPalette reliefPalette = ColorPalette::Relief();

    QImage shading(hf.getSizeX(), hf.getSizeY(), QImage::Format_ARGB32);
    for (int i = 0; i < hf.getSizeX(); i++) {
        for (int j = 0; j < hf.getSizeY(); j++) {

            Vector3 pos = hf.vertex(i, j);
            Vector2 p(pos[0], pos[1]);

            // Height shading: color is a linear interpolation of height colors
            double t = 300.0 * Math::LinearStep(hf.at(i, j), -b, 1.5 * b);
            Vector3 cz = reliefPalette.getColor(t);
            cz = 0.25 * cz;
            cz = Vector3(1 * cz[0], 0.7 * cz[1], 0.5 * cz[2]);

            // Ice shading
            Vector3 ci;
            if (ice.at(i, j) > 0.0) {
                ci = Vector3(148 / 255.0, 164 / 255.0, 164 / 255.0) * (1 - ela.value(p)) + (ela.value(p)) * Vector3(1, 1, 1);

                // blue-ish ramp
                double tb = (1 - Math::LinearStep(hf.at(i, j), 1000, 2500.0));
                Vector3 cb = Vector3(ci[0] * 0.2, ci[1], ci[2]);
                ci = (1 - tb) * ci + tb * cb;

                // thin ice: interpolate towards rock texture
                double ti = std::min(ice.at(i, j) / 2.0, 1.0);
                if (ela.value(p) < 1) ti = 1;
                cz = (1 - ti) * cz + ti * ci;
            }

            // light
            double light = dot(hf.normal(i, j), lightDir);
            light = 0.5 * (1.0 + light);
            double cs = 0.9 * light;
            cs *= cs; // cosine like

            // Normal shading: color is a combination of cool and cold colors according to the orientation
            Vector3 ambient = 0.25 * Vector3(1.0, 1.0, 1.0);
            Vector3 c1 = 0.25 * cs * cz;
            Vector3 c2 = 0.5 * cs * Vector3(1.0, 1.0, 1.0);
            Vector3 c = ambient + c1 + c2;

            shading.setPixelColor(i, j, c.toQColor().rgb());
        }
    }

    int texw = int(hf.getDomain().width() / TEXTURE_RES);
    int texh = int(hf.getDomain().height() / TEXTURE_RES);

    std::cerr << "Creating base texture of size " << texw << "x" << texh << std::endl;

    QImage tex(texw, texh, QImage::Format_RGBA8888);
    QPainter painter(&tex);
    painter.drawImage(QPoint(0, 0), shading.scaled(texw, texh, Qt::KeepAspectRatio, Qt::SmoothTransformation));
    return tex;
}


QImage FeaturesPainter::crevasses(const ScalarField2& crevasseField, int crevasseType, double scale, const ScalarField2& alphaCrevasses, QImage& texture)
{
    QImage crevasseMap(texture.size(), QImage::Format_Grayscale16);
    crevasseMap.fill(quint16(65535));
    QPainter mapPainter(&crevasseMap);

    // Stamps (we might be interested in having different crevasse patterns above/below ELA)
    unsigned int CREVASSE_STAMP_SIZE = 256;
    std::vector<std::vector<QImage> > stamps = {
        { // stamps under ELA
            QImage("data/stamps/crevasses_1_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_2_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_3_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_4_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_5_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_6_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_7_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation)
        },
        { // stamps above ELA
            QImage("data/stamps/crevasses_2_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_5_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_7_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_8_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation),
            QImage("data/stamps/crevasses_9_mask.png").scaledToWidth(CREVASSE_STAMP_SIZE, Qt::TransformationMode::SmoothTransformation)
        }
    };

    // Bombing
    for (int i = 0; i < crevasseField.getSizeX(); i++) {
        for (int j = 0; j < crevasseField.getSizeY(); j++) {

            double crevProb = crevasseField.at(i, j);
            double r = randUniform();

            if (r < crevProb)
            {
                // World space point
                Vector2 p = crevasseField.domainCoords(i, j);

                Vector2 flowDir = flowDirField.at(i, j);
                double flowDirAngle = Math::RadianToDegree(std::atan2(flowDir[1], flowDir[0]));
                bool aboveELA = accumArea.value(p) > 0.5 ? 1 : 0;

                // angle and size of crevasse
                double crevAngle, crevSize;
                switch (crevasseType) {

                    case 0: // transverse
                        crevAngle = flowDirAngle + 90;
                        crevSize = std::min(aboveELA ? 2.5 : 3.5, valleyDist.at(i, j) / 15.0);
                        crevSize /= TEXTURE_RES;
                        break;

                    case 1: // longitudinal
                        crevAngle = flowDirAngle;
                        crevSize = 3.0 / TEXTURE_RES;
                        break;

                    case 2: { // marginal
                        // are we on left or right side of glacier?
                        Vector2 gradDist = valleyDist.gradient(i, j);
                        if (Norm(gradDist) < 1e-3) continue;
                        Vector2 glacierCenterDir = Normalized(gradDist);
                        Vector2 ctrPoint = p + Norm(valleyDist.getCellSize()) * glacierCenterDir * valleyDist.value(p);
                        Vector2 ctrDown = ctrPoint + flowDirField.value(p);
                        Vector2 ctrFlowDir = ctrDown - ctrPoint;
                        double side = (p[0] - ctrPoint[0]) * ctrFlowDir[1] - (p[1] - ctrPoint[1]) * ctrFlowDir[0];
                        double sign = side < 0 ? -1 : 1;
                        crevAngle = Math::RadianToDegree(std::atan2(glacierCenterDir[1], glacierCenterDir[0])) + sign * 45.0;
                        crevSize = 3.0 / TEXTURE_RES;
                    }   break;

                    case 3: // icefall
                        crevAngle = flowDirAngle + randUniform(-20, 20) + 90.0 * (randUniform() > 0.5);
                        crevSize = 2.0 * randUniform(0.75, 1.5) / TEXTURE_RES;
                        break;

                    default:
                        crevAngle = 0;
                        crevSize = 0;
                        break;
                }

                // Texture coordinates (high res)
                int ii, jj;
                remapTexCoordinates(p, crevasseField.getDomain(), texture.size(), ii, jj);

                int stampIdx = int(randUniform() * (aboveELA ? stamps[1].size() : stamps[0].size()));

                // Draw map stamp
                mapPainter.setCompositionMode(QPainter::CompositionMode::CompositionMode_Multiply);
                mapPainter.translate(QPoint(ii, texture.height() - 1 - jj));
                mapPainter.rotate(-crevAngle);

                const QImage& stampToPlace = stamps[aboveELA ? 1 : 0][stampIdx];
                double cscale = randUniform(0.8, 3.0) * crevSize * 18.0 / CREVASSE_STAMP_SIZE * scale;
                mapPainter.drawImage(QPoint(0, 0), stampToPlace.scaled(stampToPlace.width() * cscale, stampToPlace.height() * cscale));

                mapPainter.rotate(crevAngle);
                mapPainter.translate(-QPoint(ii, texture.height() - 1 - jj));
            }
        }
    }

    // draw crevasses to texture and at the same time remove crevasses outside ice in map
    quint16* cmap = reinterpret_cast<quint16*>(crevasseMap.bits());
    for (int i = 0; i < crevasseMap.width(); i++) {
        for (int j = 0; j < crevasseMap.height(); j++) {
            double crevValue = cmap[j * crevasseMap.width() + i] / 65535.0;
            double alpha = alphaCrevasses.at(i, alphaCrevasses.getSizeY() - 1 - j);
            crevValue = alpha * crevValue + (1 - alpha) * 1.0;
            cmap[j * crevasseMap.width() + i] = quint16(65535 * crevValue);

            Vector3 texColor = Vector3::fromQColor(texture.pixel(i, j)) * crevValue;
            texture.setPixelColor(i, j, texColor.toQColor());
        }
    }

    return crevasseMap;
}


QImage FeaturesPainter::rimayes(const ScalarField2& rimayeField, QImage& texture)
{
    const QPoint next[8] = { QPoint(1, 0), QPoint(1, 1), QPoint(0, 1), QPoint(-1, 1), QPoint(-1, 0), QPoint(-1, -1), QPoint(0, -1), QPoint(1, -1) };

    QImage rimayesMap(texture.size(), QImage::Format_Grayscale16);
    rimayesMap.fill(quint16(65535));
    QPainter mapPainter(&rimayesMap);
    QPainter painter(&texture);

    double penSize = 3.0 / TEXTURE_RES;
    QPen rimayePenBlack(QColor(48, 48, 48), penSize, Qt::SolidLine, Qt::SquareCap, Qt::MiterJoin);

    // Drawing rimayes
    int rimayeCount = 0;
    std::vector<bool> visited(rimayeField.getSizeX() * rimayeField.getSizeY(), false);
    for (int i = 0; i < rimayeField.getSizeX(); i++) {
        for (int j = 0; j < rimayeField.getSizeY(); j++) {
            if (rimayeField.at(i, j) > 0.5 && !visited[rimayeField.vertexIndex(i, j)])
            {
                // count neighbors
                int numNeighs = 0;
                for (int n = 0; n < 8; n++) {
                    QPoint pn = QPoint(i, j) + next[n];
                    if (rimayeField.validIndex(pn.x(), pn.y()) && rimayeField.at(pn.x(), pn.y()) > 0.5
                        && !visited[rimayeField.vertexIndex(pn.x(), pn.y())])
                    {
                        numNeighs++;
                    }
                }

                // if more than one neighbor, this is an intermediate rimaye node and we will
                // visit it eventually starting from another point
                // this simplifies the logic below
                if (numNeighs > 1) continue;

                // starting at the node, visit neighbors until reaching the other rimaye endpoint
                std::vector<QPoint> rimayeNodes;
                QPoint pcurr(i, j);
                bool foundNext = true;
                while (foundNext) {
                    rimayeNodes.push_back(pcurr);
                    visited[rimayeField.vertexIndex(pcurr.x(), pcurr.y())] = true;
                    foundNext = false;
                    for (int n = 0; n < 8; n++) {
                        QPoint pnext = pcurr + next[n];
                        if (rimayeField.validIndex(pnext.x(), pnext.y()) && rimayeField.at(pnext.x(), pnext.y()) > 0.5
                            && !visited[rimayeField.vertexIndex(pnext.x(), pnext.y())]) {
                            pcurr = pnext;
                            foundNext = true;
                            break;
                        }
                    }
                }

                if (rimayeNodes.size() < 2) continue;

                // draw rimaye
                bool first = true;
                QPainterPath path;
                for (QPoint pn : rimayeNodes) {
                    // World space point
                    Vector2 p = rimayeField.domainCoords(pn.x(), pn.y()) + 0.25 * randUniform() * rimayeField.getCellSize();

                    // Texture coordinates (high res)
                    int ii, jj;
                    remapTexCoordinates(p, rimayeField.getDomain(), texture.size(), ii, jj);
                    if (first) {
                        path.moveTo(ii, texture.height() - 1 - jj);
                        first = false;
                    }
                    else {
                        path.lineTo(ii, texture.height() - 1 - jj);
                    }
                }

                painter.setPen(rimayePenBlack);
                painter.drawPath(path);

                mapPainter.setPen(rimayePenBlack);
                mapPainter.drawPath(path);

                rimayeCount++;
            }
        }
    }

    return rimayesMap;
}


QImage FeaturesPainter::ogives(const ScalarField2& ogiveDistSrc, const ScalarField2& ogiveDistCtr, const QImage& ogiveId, double ogiveLenScale, double ogiveFreq,
    QImage& texture, ScalarField2& ogiveStrength)
{
    QImage ogiveMap(texture.size(), QImage::Format_Grayscale16);
    ogiveMap.fill(32767);
    quint16* ogiveMapData = (quint16*)ogiveMap.bits();

    std::vector<double> ogiveLength(255, -1);
    std::vector<double> ogiveScale(255, -1);

    for (int i = 0; i < texture.width(); i++) {
        for (int j = 0; j < texture.height(); j++) {

            // world coords
            double di = i / double(texture.width());
            double dj = j / double(texture.height());
            double x = di * ogiveDistSrc.getDomain().width();
            double y = dj * ogiveDistSrc.getDomain().height();
            Vector2 p = Vector2(x + ogiveDistSrc.getDomain().getMin()[0],
                ogiveDistSrc.getDomain().getMax()[1] - y);

            // small texture coords
            int ii = int(di * ogiveId.width());
            int jj = int(dj * ogiveId.height());
            int oid = qRed(ogiveId.pixel(ii, ogiveId.height() - 1 - jj));
            if (oid > 0) {

                if (ogiveLength[oid] < 0) {
                    ogiveLength[oid] = randUniform(2000, 5000) * ogiveLenScale;
                    ogiveScale[oid] = randUniform(0.3, 1.0);
                }

                // create waves pattern based on dist
                double distSrc = ogiveDistSrc.value(p);
                double distCtr = ogiveDistCtr.value(p);
                double dist = distSrc + distCtr;
                double wave = std::cos(dist*ogiveFreq / 15.0);
                double damp = std::max(0.0, 1.0 - std::max(0.0, dist - 0.2 * ogiveLength[oid]) / ogiveLength[oid]);

                double hOff = ogiveScale[oid] * wave * damp;
                int    cOff = int(255 * 0.05 * damp * wave);

                QColor c = texture.pixelColor(i, j);
                c = QColor(std::min(c.red() + cOff, 255),
                    std::min(c.green() + cOff, 255),
                    std::min(c.blue() + cOff, 255),
                    c.alpha());
                texture.setPixelColor(i, j, c);

                ogiveMapData[i + j * texture.width()] = quint16(std::min(std::max(0.0, (0.5 * hOff + 0.5) * 65535), 65535.0));

                ogiveStrength(i, ogiveStrength.getSizeY() - 1 - j) = damp;
            }
        }
    }

    return ogiveMap;
}


QImage FeaturesPainter::moraines(const Vector3& moraineColor, double moraineScale, QImage& texture)
{
    const unsigned int minPathNodes = 100;

    // upscale map using NN to avoid one-pixel segmentations
    QImage glaciersIdUp = glacierIds.scaled(QSize(2 * glacierIds.width(), 2 * glacierIds.height()), Qt::KeepAspectRatio, Qt::FastTransformation);

    // 1: identify all of the nodes of a moraine based on segmentation map
    std::map<std::pair<int, int>, std::vector<Vector2> > moraineNodes;
    for (int i = 1; i < glaciersIdUp.width() - 1; i++) {
        for (int j = 1; j < glaciersIdUp.height() - 1; j++) {
            // fetch ids
            int glacierId = qBlue(glaciersIdUp.pixel(i, glaciersIdUp.height() - 1 - j));
            int subflowId = qGreen(glaciersIdUp.pixel(i, glaciersIdUp.height() - 1 - j));
            int gid = glacierId * 256 + subflowId;
            if (glacierId == 0 || subflowId == 0) continue;

            // check neighbors
            for (int di = 0; di <= 1; di++) {
                for (int dj = 0; dj <= 1; dj++) {
                    if (di * dj > 0) continue;
                    int neighGlacier = qBlue(glaciersIdUp.pixel(i + di, glaciersIdUp.height() - 1 - j - dj));
                    int neighSubflow = qGreen(glaciersIdUp.pixel(i + di, glaciersIdUp.height() - 1 - j - dj));
                    int neighId = neighGlacier * 256 + neighSubflow;

                    // two glacier flows touching
                    if (neighGlacier == glacierId && neighSubflow != subflowId && neighSubflow > 0) {
                        std::pair<int, int> idPair(std::min(neighId, gid), std::max(neighId, gid));
                        if (moraineNodes.find(idPair) == moraineNodes.end()) {
                            moraineNodes[idPair] = std::vector<Vector2>();
                        }
                        moraineNodes[idPair].push_back(Vector2(i + 0.5 + 0.5 * di, j + 0.5 + 0.5 * dj));
                    }
                }
            }
        }
    }

    // 2: create moraines paths
    std::vector<std::vector<Vector2> > morainePaths;
    std::vector<std::vector<Vector2> > partialPaths;
    for (auto it = moraineNodes.begin(); it != moraineNodes.end(); it++) {

        // find all node neighbors pairs
        std::vector<Vector2> nodes = it->second;
        std::vector<std::vector<int> > nodeNeighs(nodes.size());
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                if (SquaredNorm(nodes[i] - nodes[j]) < 1.01) {
                    nodeNeighs[i].push_back(j);
                    nodeNeighs[j].push_back(i);
                }
            }
        }

        // starting from a node with only 1 neighbor, propagate until end
        std::vector<bool> nodeSeen(nodes.size(), false);
        for (int i = 0; i < nodes.size(); i++) {
            if (nodeSeen[i]) continue;
            if (nodeNeighs[i].size() == 1) {

                // build path
                std::vector<Vector2> path;
                int prev = -1;
                int curr = i;
                while (curr >= 0) {
                    path.push_back(nodes[curr]);
                    nodeSeen[curr] = true;
                    int next = -1;
                    for (int j : nodeNeighs[curr]) {
                        if (j != prev && !nodeSeen[j]) {
                            next = j;
                            break;
                        }
                    }
                    prev = curr;
                    curr = next;
                }

                if (path.size() > minPathNodes) {
                    morainePaths.push_back(path);
                }
                else if (path.size() > 4) {
                    partialPaths.push_back(path);
                }
            }
        }
    }

    // complete paths with partial segments
    int addedPaths = 1;
    std::vector<bool> partialPathAdded(partialPaths.size(), false);
    while (addedPaths > 0) {
        addedPaths = 0;
        for (int pi = 0; pi < partialPaths.size(); pi++) {
            if (partialPathAdded[pi]) continue;
            const std::vector<Vector2>& pp = partialPaths[pi];
            Vector2 partialFront = pp.front();
            Vector2 partialBack = pp.back();
            for (int mj = 0; mj < morainePaths.size(); mj++) {
                std::vector<Vector2>& path = morainePaths[mj];
                Vector2 pathFront = path.front();
                Vector2 pathBack = path.back();

                const int maxDist = 10;//3;
                if (SquaredNorm(partialFront - pathFront) < maxDist) {
                    for (int i = 0; i < pp.size(); i++) {
                        path.insert(path.begin(), pp[i]);
                    }
                    addedPaths++;
                    partialPathAdded[pi] = true;
                    break;
                }
                if (SquaredNorm(partialFront - pathBack) < maxDist) {
                    for (int i = 0; i < pp.size(); i++) {
                        path.insert(path.end(), pp[i]);
                    }
                    addedPaths++;
                    partialPathAdded[pi] = true;
                    break;
                }
                if (SquaredNorm(partialBack - pathFront) < maxDist) {
                    for (int i = int(pp.size()) - 1; i >= 0; i--) {
                        path.insert(path.begin(), pp[i]);
                    }
                    addedPaths++;
                    partialPathAdded[pi] = true;
                    break;
                }
                if (SquaredNorm(partialBack - pathBack) < maxDist) {
                    for (int i = int(pp.size()) - 1; i >= 0; i--) {
                        path.insert(path.end(), pp[i]);
                    }
                    addedPaths++;
                    partialPathAdded[pi] = true;
                    break;
                }
            }
        }
    }

    // smoothing of positions
    const int SMOOTH_ITERS = 10;
    for (int s = 0; s < SMOOTH_ITERS; s++) {
        for (int m = 0; m < morainePaths.size(); m++) {
            std::vector<Vector2> newPositions = morainePaths[m];
            for (int i = 1; i < newPositions.size() - 1; i++) {
                newPositions[i] = 0.5 * (morainePaths[m][i - 1] + morainePaths[m][i + 1]);
            }
            morainePaths[m] = newPositions;
        }
    }

    // 3: draw moraines
    const int moraineMapDownscale = 2;
    QImage moraineMap(QSize(texture.width() / moraineMapDownscale, texture.height() / moraineMapDownscale), QImage::Format_Grayscale16);
    moraineMap.fill(0);
    QPen morainePen(QColor(255, 255, 255), 1, Qt::SolidLine, Qt::PenCapStyle::FlatCap, Qt::PenJoinStyle::MiterJoin);
    QPainter painterMap(&moraineMap);

    for (int m = 0; m < morainePaths.size(); m++) {
        const auto moraine = morainePaths[m];

        // go to tex coords
        std::vector<QPoint> path;
        for (const Vector2& p : moraine) {
            double di = p[0] / double(glaciersIdUp.width());
            double dj = p[1] / double(glaciersIdUp.height());
            int ii = int(di * texture.width() / moraineMapDownscale + 0.5);
            int jj = int(dj * texture.height() / moraineMapDownscale + 0.5);
            path.push_back(QPoint(ii, jj));
        }

        QPainterPath ppath;
        ppath.moveTo(path[0]);
        for (int p = 0; p < path.size(); p++) {
            ppath.lineTo(path[p]);
        }
        morainePen.setWidth(randUniform(5.0, 50.0) / TEXTURE_RES * moraineScale);
        painterMap.setPen(morainePen);
        painterMap.drawPath(ppath);
    }

    // 4: round shape and smooth fading
    ScalarField2 morainesField = ScalarField2(terrainBox, moraineMap, 0.0, 1.0, true);
    for (int blurPass = 0; blurPass < 4; blurPass++)
        morainesField.gaussianBlur();
    //morainesField.Normalize();
    morainesField = morainesField.setResolution(texture.width(), texture.height());

    for (int i = 0; i < morainesField.getSizeX(); i++) {
        for (int j = 0; j < morainesField.getSizeY(); j++) {

            // world coords
            double di = i / double(texture.width());
            double dj = j / double(texture.height());
            double x = di * terrainBox.width();
            double y = dj * terrainBox.height();
            Vector2 p = Vector2(x + terrainBox.getMin()[0], terrainBox.getMax()[1] - y);

            // fetch ids
            //int ii = int(di * glaciersIdUp.width());
            //int jj = int(dj * glaciersIdUp.height());
            //int glacierId = qBlue(glaciersIdUp.pixel(ii, glaciersIdUp.height() - 1 - jj));
            //int subflowId = qGreen(glaciersIdUp.pixel(ii, glaciersIdUp.height() - 1 - jj));

            double moraine = morainesField.at(i, j);
            moraine = std::sqrt(moraine);
            // remove moraines above ELA
            moraine *= 1.0 - accumArea.value(p);
            // fade out moraines on very steep areas
            moraine *= 1.0 - std::min(1.0, std::max(0.0, slopeField.value(p) - 0.2) / 0.3);
            //moraine *= subflowId > 0;
            morainesField(i, j) = moraine;

            // final color
            Vector3 c = Vector3::fromQColor(texture.pixelColor(i, j));
            Vector3 cfinal = moraine * moraineColor + (1 - moraine) * c;
            texture.setPixelColor(i, j, cfinal.toQColor());
        }
    }

    QImage resultMap(QSize(texture.width(), texture.height()), QImage::Format_Grayscale16);
    quint16* resData = (quint16*)resultMap.bits();
    for (int i = 0; i < resultMap.width(); i++) {
        for (int j = 0; j < resultMap.height(); j++) {
            resData[i + j * resultMap.width()] = quint16(65535 * morainesField.at(i, j));
        }
    }

    return resultMap;
}


ScalarField2 FeaturesPainter::seracs(const ScalarField2& seracs, double displacement)
{
    ScalarField2 seracIce(iceHi);

    SimplexNoise simplex;

    int nx = iceLow.getSizeX();
    int ny = iceLow.getSizeY();
    int hx = iceHi.getSizeX();
    int hy = iceHi.getSizeY();

    ScalarField2 downIce2 = iceLow.setResolution(nx / 2, ny / 2);
    ScalarField2 downIce4 = iceLow.setResolution(nx / 4, ny / 4);
    ScalarField2 downIce8 = iceLow.setResolution(nx / 8, ny / 8);

    for (int i = 0; i < hx; i++) {
        for (int j = 0; j < hy; j++) {
            // low coordinates
            double x = double(i) / hx;
            double y = double(j) / hy;
            int ii = int(x * nx);
            int jj = int(y * ny);

            double hnoise = -2.5 + 5.0 * Math::Ridge(simplex.at(30.0 * i / 250.0, 30.0 * j / 250.0), 0.5);
            double loIce2 = downIce2.at(int(x * nx / 2), int(y * ny / 2)) - iceLow.at(ii, jj);
            double loIce4 = downIce4.at(int(x * nx / 4), int(y * ny / 4)) - iceLow.at(ii, jj);
            double loIce8 = downIce8.at(int(x * nx / 8), int(y * ny / 8)) - iceLow.at(ii, jj);

            double iceDisp = (1.5 * hnoise + 1.0 * loIce8 + 0.5 * loIce4 + 0.25 * loIce2);
            double iceL = std::min(iceLow.at(ii,jj) + iceDisp*displacement, iceHi.at(i, j));

            double t = seracs.at(ii, jj);
            double ice = iceL * t + (1 - t) * iceHi.at(i, j);
            seracIce(i, j) = ice;
        }
    }

    return seracIce;
}
