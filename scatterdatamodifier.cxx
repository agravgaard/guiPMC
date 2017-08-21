/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "scatterdatamodifier.h"
#include <QtDataVisualization/qscatterdataproxy.h>
#include <QtDataVisualization/qvalue3daxis.h>
#include <QtDataVisualization/q3dscene.h>
#include <QtDataVisualization/q3dcamera.h>
#include <QtDataVisualization/qscatter3dseries.h>
#include <QtDataVisualization/q3dtheme.h>
#include <QtCore/qmath.h>
#include <QtWidgets/QComboBox>

using namespace QtDataVisualization;

//#define RANDOM_SCATTER // Uncomment this to switch to random scatter

const int numberOfItems = 3600;
const float curveDivider = 3.0f;
const int lowerNumberOfItems = 900;
const float lowerCurveDivider = 0.75f;

ScatterDataModifier::ScatterDataModifier(Q3DScatter *scatter)
	: m_graph(scatter),
	m_style(QAbstract3DSeries::MeshSphere),
	m_smooth(true),
	m_itemCount(lowerNumberOfItems),
	m_curveDivider(lowerCurveDivider)
{
	m_graph->activeTheme()->setType(Q3DTheme::ThemeQt);
	m_graph->setShadowQuality(QAbstract3DGraph::ShadowQualityNone);
	m_graph->scene()->activeCamera()->setCameraPreset(Q3DCamera::CameraPresetFront);
	// m_graph->setOptimizationHints(QAbstract3DGraph::OptimizationStatic); // Optimize for performance for ref: ( http://blog.qt.io/blog/2015/02/10/qt-weekly-25-optimization-hints-on-qt-data-visualization/ )
	// ^causes allocation size assertion exception in qVector.cpp at draw-time
	QScatterDataProxy *proxy = new QScatterDataProxy;
	QScatter3DSeries *series = new QScatter3DSeries(proxy);
	series->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
	series->setMeshSmooth(m_smooth);
	m_graph->addSeries(series);

	m_pos = NULL;
	m_dir = NULL;
	m_weight = NULL;
	m_energy = NULL;
#ifdef QT_TEST
	addData();
#endif
}

ScatterDataModifier::~ScatterDataModifier()
{
	delete m_graph; // is this recursive to all members?
	delete[] m_pos;
	delete[] m_dir;
	delete[] m_weight;
	delete[] m_energy;
}

void ScatterDataModifier::setItemCount(size_t itemcount){
	m_itemCount = itemcount;
}

void ScatterDataModifier::setPosData(cl_float3* pos){
	if (!m_pos == NULL)
		delete[] m_pos;

	m_pos = pos;
}
void ScatterDataModifier::setDirectionData(cl_float3* dir){
	if (!m_dir == NULL)
		delete[] m_dir;

	m_dir = dir;
}
void ScatterDataModifier::setWeightData(cl_float* weight){
	if (!m_weight == NULL)
		delete[] m_weight;

	m_weight = weight;
}
void ScatterDataModifier::setEnergyData(cl_float* T){
	if (!m_energy == NULL)
		delete[] m_energy;

	m_energy = T;
}

int get_idx(const float local_weight, const float color_step, const float min_weight, const size_t N_color_steps){
	int i = 0;
	while (i < N_color_steps){
		if ((min_weight + i*color_step) < local_weight && local_weight < (min_weight + (i + 1)*color_step))
			return i;
		i++;
	}
	return N_color_steps - 1;
}

void ScatterDataModifier::addData()
{
	// Configure the axes according to the data
	m_graph->axisX()->setTitle("X");
	m_graph->axisY()->setTitle("Y");
	m_graph->axisZ()->setTitle("Z");

#ifdef QT_TEST
#ifdef RANDOM_SCATTER
	for (int i = 0; i < m_itemCount; i++) {
		ptrToDataArray->setPosition(randVector());
		ptrToDataArray++;
	}
#else
	float limit = qSqrt(m_itemCount) / 2.0f;
	for (float i = -limit; i < limit; i++) {
		for (float j = -limit; j < limit; j++) {
			ptrToDataArray->setPosition(QVector3D(i + 0.5f,
				qCos(qDegreesToRadians((i * j) / m_curveDivider)),
				j + 0.5f));
			ptrToDataArray++;
		}
	}
#endif
#else
	float min_weigth = 99999.9;
	float max_weigth = -1.0;
	int limit = m_itemCount;

	// calculate color weighting scale
	for (int i = 0; i < limit; i++) {
		const float local_weight = m_energy[i];
		if (local_weight > max_weigth)
			max_weigth = local_weight;
		if (local_weight < min_weigth)
			min_weigth = local_weight;
	}
	// make N color steps:
	std::vector<QColor> color_gradient;
	color_gradient.push_back(QColor(QRgb(0x0000ff)));
	color_gradient.push_back(QColor(QRgb(0x0033ff)));
	color_gradient.push_back(QColor(QRgb(0x0077ff)));
	color_gradient.push_back(QColor(QRgb(0x00ff77)));
	color_gradient.push_back(QColor(QRgb(0x00ff33)));
	color_gradient.push_back(QColor(QRgb(0x00ff00)));
	color_gradient.push_back(QColor(QRgb(0x33ff00)));
	color_gradient.push_back(QColor(QRgb(0x77ff00)));
	color_gradient.push_back(QColor(QRgb(0xff7700)));
	color_gradient.push_back(QColor(QRgb(0xff3300)));
	color_gradient.push_back(QColor(QRgb(0xff0000)));

	const size_t N_color_steps = color_gradient.size();
	const float color_step = (max_weigth - min_weigth) / float(N_color_steps);

	std::vector<QScatterDataArray *> dataArray;
	std::vector<QScatterDataItem *> ptrToDataArray;
	for (int idx = 0; idx < N_color_steps; idx++){
		dataArray.push_back(new QScatterDataArray);
		dataArray.at(idx)->resize(m_itemCount);
		ptrToDataArray.push_back(&dataArray.at(idx)->first());
	}

	for (int i = 0; i < limit; i++) {
		const float local_weigth = m_weight[i] * 50.0f;
		const int idx = get_idx(m_energy[i], color_step, min_weigth, N_color_steps);
		ptrToDataArray.at(idx)->setPosition(QVector3D(
			m_pos[i].x + m_dir[i].x * local_weigth, // x pos (floor 1st axis)
			m_pos[i].y + m_dir[i].y * local_weigth, // y pos (up in viewer)
			m_pos[i].z + m_dir[i].z * local_weigth)); // z pos (floor 2nd axis)
		ptrToDataArray.at(idx)++;
	}
#endif

	m_graph->seriesList().at(0)->dataProxy()->resetArray(dataArray.at(0));
	m_graph->seriesList().at(0)->setBaseColor(color_gradient.at(0));

	for (int idx = 1; idx < N_color_steps; idx++){
		QScatterDataProxy* new_dataproxy = new QScatterDataProxy(m_graph->seriesList().at(0)->dataProxy()->parent());
		new_dataproxy->resetArray(dataArray.at(idx));
		QScatter3DSeries* new_series = new QScatter3DSeries(new_dataproxy, m_graph->seriesList().at(0)->parent());
		new_series->setBaseColor(color_gradient.at(idx));
		new_series->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
		new_series->setMeshSmooth(m_smooth);
		m_graph->addSeries(new_series);
	}
}

void ScatterDataModifier::callAddData()
{
	addData();
}

void ScatterDataModifier::changeStyle(int style)
{
	QComboBox *comboBox = qobject_cast<QComboBox *>(sender());
	if (comboBox) {
		m_style = QAbstract3DSeries::Mesh(comboBox->itemData(style).toInt());
		for (int i = 0; i < m_graph->seriesList().size(); i++)
			m_graph->seriesList().at(i)->setMesh(m_style);
	}
}

void ScatterDataModifier::setSmoothDots(int smooth)
{
	m_smooth = bool(smooth);
	for (int i = 0; i < m_graph->seriesList().size(); i++){
		QScatter3DSeries *series = m_graph->seriesList().at(i);
		series->setMeshSmooth(m_smooth);
	}
}

void ScatterDataModifier::changeTheme(int theme)
{
	Q3DTheme *currentTheme = m_graph->activeTheme();
	currentTheme->setType(Q3DTheme::Theme(theme));
	emit backgroundEnabledChanged(currentTheme->isBackgroundEnabled());
	emit gridEnabledChanged(currentTheme->isGridEnabled());
	emit fontChanged(currentTheme->font());
}

void ScatterDataModifier::changePresetCamera()
{
	static int preset = Q3DCamera::CameraPresetFrontLow;

	m_graph->scene()->activeCamera()->setCameraPreset((Q3DCamera::CameraPreset)preset);

	if (++preset > Q3DCamera::CameraPresetDirectlyBelow)
		preset = Q3DCamera::CameraPresetFrontLow;
}

void ScatterDataModifier::changeLabelStyle()
{
	m_graph->activeTheme()->setLabelBackgroundEnabled(!m_graph->activeTheme()->isLabelBackgroundEnabled());
}

void ScatterDataModifier::shadowQualityUpdatedByVisual(QAbstract3DGraph::ShadowQuality sq)
{
	int quality = int(sq);
	emit shadowQualityChanged(quality); // connected to a checkbox in main.cpp
}

void ScatterDataModifier::changeShadowQuality(int quality)
{
	QAbstract3DGraph::ShadowQuality sq = QAbstract3DGraph::ShadowQuality(quality);
	m_graph->setShadowQuality(sq);
}

void ScatterDataModifier::setBackgroundEnabled(int enabled)
{
	m_graph->activeTheme()->setBackgroundEnabled((bool)enabled);
}

void ScatterDataModifier::setGridEnabled(int enabled)
{
	m_graph->activeTheme()->setGridEnabled((bool)enabled);
}

#ifdef QT_TEST
void ScatterDataModifier::toggleItemCount()
{
	if (m_itemCount == numberOfItems) {
		m_itemCount = lowerNumberOfItems;
		m_curveDivider = lowerCurveDivider;
	}
	else {
		m_itemCount = numberOfItems;
		m_curveDivider = curveDivider;
	}
	m_graph->seriesList().at(0)->dataProxy()->resetArray(0);
	addData();
}

QVector3D ScatterDataModifier::randVector()
{
	return QVector3D(
		(float)(rand() % 100) / 2.0f - (float)(rand() % 100) / 2.0f,
		(float)(rand() % 100) / 100.0f - (float)(rand() % 100) / 100.0f,
		(float)(rand() % 100) / 2.0f - (float)(rand() % 100) / 2.0f);
}
#endif