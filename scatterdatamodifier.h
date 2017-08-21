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

#ifndef SCATTERDATAMODIFIER_H
#define SCATTERDATAMODIFIER_H

#include <QtDataVisualization/q3dscatter.h>
#include <QtDataVisualization/qabstract3dseries.h>
#include "cl.hpp"

using namespace QtDataVisualization;

class ScatterDataModifier : public QObject
{
	Q_OBJECT
public:
	explicit ScatterDataModifier(Q3DScatter *scatter);
	~ScatterDataModifier();

	void setItemCount(size_t itemcount);
	void setPosData(cl_float3* pos);
	void setDirectionData(cl_float3* dir);
	void setWeightData(cl_float* weight);
	void setEnergyData(cl_float* T);
	void addData();
	void callAddData();
	void changeStyle();
	void changePresetCamera();
	void changeLabelStyle();
	void setBackgroundEnabled(int enabled);
	void setGridEnabled(int enabled);
	void setSmoothDots(int smooth);
	void toggleItemCount();
	void start();

	public Q_SLOTS:
	void changeStyle(int style);
	void changeTheme(int theme);
	void changeShadowQuality(int quality);
	void shadowQualityUpdatedByVisual(QAbstract3DGraph::ShadowQuality shadowQuality);

Q_SIGNALS:
	void backgroundEnabledChanged(bool enabled);
	void gridEnabledChanged(bool enabled);
	void shadowQualityChanged(int quality);
	void fontChanged(QFont font);

private:
	QVector3D randVector();
	Q3DScatter *m_graph;
	QAbstract3DSeries::Mesh m_style;
	bool m_smooth;
	int m_itemCount;
	cl_float3* m_pos;
	cl_float3* m_dir;
	cl_float* m_weight;
	cl_float* m_energy;
	float m_curveDivider;
};

#endif