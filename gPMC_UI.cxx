/*
*
* file  gPMC_UI.cpp
* UI program for creating, inspecting and visualizing goPMC from commands.
*
* author Andreas Gravgaard Andersen
*
* last update on 26/7/2017
*
*/
#include "gPMC_UI.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

// ITK stuff
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"

#include <QFileInfo>
#include <QFileDialog>
#include <QProcess>

#include "gPMC_dcm_tools.hxx"

gPMC_UI::gPMC_UI(QWidget *parent, Qt::WindowFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);
	connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Move()), this, SLOT(SLT_UpdateSplit1())); //added
	connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Move()), this, SLOT(SLT_UpdateSplit2())); //added
	connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Move()), this, SLOT(SLT_UpdateSplit3())); //added

	connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Pressed_Left()), this, SLOT(SLT_MousePressedLeft1())); //added
	connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Pressed_Left()), this, SLOT(SLT_MousePressedLeft2())); //added
	connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Pressed_Left()), this, SLOT(SLT_MousePressedLeft3())); //added
	connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Pressed_Right()), this, SLOT(SLT_MousePressedRight1())); //added
	connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Pressed_Right()), this, SLOT(SLT_MousePressedRight2())); //added
	connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Pressed_Right()), this, SLOT(SLT_MousePressedRight3())); //added

	connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Released_Left()), this, SLOT(SLT_MouseReleasedLeft1())); //added
	connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Released_Left()), this, SLOT(SLT_MouseReleasedLeft2())); //added
	connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Released_Left()), this, SLOT(SLT_MouseReleasedLeft3())); //added
	connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Released_Right()), this, SLOT(SLT_MouseReleasedRight1())); //added
	connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Released_Right()), this, SLOT(SLT_MouseReleasedRight2())); //added
	connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Released_Right()), this, SLOT(SLT_MouseReleasedRight3())); //added

	connect(ui.labelOverlapWnd1, SIGNAL(FocusOut()), this, SLOT(SLT_CancelMouseAction())); //added
	connect(ui.labelOverlapWnd2, SIGNAL(FocusOut()), this, SLOT(SLT_CancelMouseAction())); //added
	connect(ui.labelOverlapWnd3, SIGNAL(FocusOut()), this, SLOT(SLT_CancelMouseAction())); //added

	connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Wheel()), this, SLOT(SLT_MouseWheelUpdate1())); //added
	connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Wheel()), this, SLOT(SLT_MouseWheelUpdate2())); //added
	connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Wheel()), this, SLOT(SLT_MouseWheelUpdate3())); //added

	SLT_CancelMouseAction();
	m_enViewArrange = AXIAL_FRONTAL_SAGITTAL;

	m_ptTmpOriginalDataOffset = QPoint(0, 0);

	m_iTmpOriginalW = 0;
	m_iTmpOriginalL = 0;
}

gPMC_UI::~gPMC_UI()
{
	for (int i = 0; i < 3; i++)
	{
		m_YKImgFixed[i].ReleaseBuffer();
		m_YKImgMoving[i].ReleaseBuffer();
		m_YKDisp[i].ReleaseBuffer();
		m_DoseImgFixed[i].ReleaseBuffer();
		m_DoseImgMoving[i].ReleaseBuffer();
		m_AGDisp_Overlay[i].ReleaseBuffer();
	}
}

void gPMC_UI::SLT_CrntPosGo()
{
	if (!m_spFixed)
		return;

	// DICOMN position, mm
	double curDCMPosX = ui.lineEditCurPosX->text().toDouble();
	double curDCMPosY = ui.lineEditCurPosY->text().toDouble();
	double curDCMPosZ = ui.lineEditCurPosZ->text().toDouble();

	UShortImageType::PointType imgOrigin = m_spFixed->GetOrigin();
	UShortImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();

	int iSliderPosIdxZ = qRound((curDCMPosZ - imgOrigin[2]) / (double)imgSpacing[2]);
	int iSliderPosIdxY = qRound((curDCMPosY - imgOrigin[1]) / (double)imgSpacing[1]);
	int iSliderPosIdxX = qRound((curDCMPosX - imgOrigin[0]) / (double)imgSpacing[0]);

	ui.sliderPosDisp1->setValue(iSliderPosIdxZ);
	ui.sliderPosDisp2->setValue(iSliderPosIdxY);
	ui.sliderPosDisp3->setValue(iSliderPosIdxX);
}

void gPMC_UI::SLT_DrawImageWhenSliceChange()
{
	if (!m_spFixed)
		return;

	int sliderPosIdxZ, sliderPosIdxY, sliderPosIdxX;

	switch (m_enViewArrange)
	{
	case AXIAL_FRONTAL_SAGITTAL:
		sliderPosIdxZ = ui.sliderPosDisp1->value();   // Z corresponds to axial, Y to frontal, X to sagittal
		sliderPosIdxY = ui.sliderPosDisp2->value();
		sliderPosIdxX = ui.sliderPosDisp3->value();
		break;
	case FRONTAL_SAGITTAL_AXIAL:
		sliderPosIdxY = ui.sliderPosDisp1->value();
		sliderPosIdxX = ui.sliderPosDisp2->value();
		sliderPosIdxZ = ui.sliderPosDisp3->value();

		break;
	case SAGITTAL_AXIAL_FRONTAL:
		sliderPosIdxX = ui.sliderPosDisp1->value();
		sliderPosIdxZ = ui.sliderPosDisp2->value();
		sliderPosIdxY = ui.sliderPosDisp3->value();
		break;
	default:
		sliderPosIdxZ = ui.sliderPosDisp1->value();   // Z corresponds to axial, Y to frontal, X to sagittal
		sliderPosIdxY = ui.sliderPosDisp2->value();
		sliderPosIdxX = ui.sliderPosDisp3->value();
		break;
	}

	UShortImageType::SizeType imgSize = m_spFixed->GetRequestedRegion().GetSize(); //1016x1016 x z
	UShortImageType::PointType imgOrigin = m_spFixed->GetOrigin();
	UShortImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();

	double curPhysPos[3];
	curPhysPos[0] = imgOrigin[2] + sliderPosIdxZ*imgSpacing[2]; //Z in default setting
	curPhysPos[1] = imgOrigin[1] + sliderPosIdxY*imgSpacing[1]; //Y
	curPhysPos[2] = imgOrigin[0] + sliderPosIdxX*imgSpacing[0]; //Z

	int refIdx = 3 - m_enViewArrange;

	if (ui.checkBoxDrawCrosshair->isChecked())
	{
		m_YKDisp[refIdx % 3].m_bDrawCrosshair = true;
		m_YKDisp[(refIdx + 1) % 3].m_bDrawCrosshair = true;
		m_YKDisp[(refIdx + 2) % 3].m_bDrawCrosshair = true;

		m_YKDisp[refIdx % 3].m_ptCrosshair.setX(sliderPosIdxX); //axial
		m_YKDisp[refIdx % 3].m_ptCrosshair.setY(sliderPosIdxY);

		m_YKDisp[(refIdx + 1) % 3].m_ptCrosshair.setX(sliderPosIdxX); //Frontal
		m_YKDisp[(refIdx + 1) % 3].m_ptCrosshair.setY(imgSize[2] - sliderPosIdxZ - 1);

		m_YKDisp[(refIdx + 2) % 3].m_ptCrosshair.setX(sliderPosIdxY); //Sagittal
		m_YKDisp[(refIdx + 2) % 3].m_ptCrosshair.setY(imgSize[2] - sliderPosIdxZ - 1);

		m_YKImgFixed[0].m_bDrawCrosshair = true;
		m_YKImgFixed[1].m_bDrawCrosshair = true;
		m_YKImgFixed[2].m_bDrawCrosshair = true;

		m_YKImgFixed[0].m_ptCrosshair.setX(sliderPosIdxX); //sagittal slider
		m_YKImgFixed[0].m_ptCrosshair.setY(sliderPosIdxY);

		m_YKImgFixed[1].m_ptCrosshair.setX(sliderPosIdxX); //sagittal slider
		m_YKImgFixed[1].m_ptCrosshair.setY(imgSize[2] - sliderPosIdxZ - 1);

		m_YKImgFixed[2].m_ptCrosshair.setX(sliderPosIdxY); //sagittal slider
		m_YKImgFixed[2].m_ptCrosshair.setY(imgSize[2] - sliderPosIdxZ - 1);
	}
	else
	{
		m_YKDisp[0].m_bDrawCrosshair = false;
		m_YKDisp[1].m_bDrawCrosshair = false;
		m_YKDisp[2].m_bDrawCrosshair = false;

		m_YKImgFixed[0].m_bDrawCrosshair = false;
		m_YKImgFixed[1].m_bDrawCrosshair = false;
		m_YKImgFixed[2].m_bDrawCrosshair = false;
	}

	if (m_spMoving)
	{
		this->Draw2DFrom3DDouble(m_spFixed, m_spMoving, PLANE_AXIAL, curPhysPos[0], m_YKImgFixed[0], m_YKImgMoving[0]);
		this->Draw2DFrom3DDouble(m_spFixed, m_spMoving, PLANE_FRONTAL, curPhysPos[1], m_YKImgFixed[1], m_YKImgMoving[1]);
		this->Draw2DFrom3DDouble(m_spFixed, m_spMoving, PLANE_SAGITTAL, curPhysPos[2], m_YKImgFixed[2], m_YKImgMoving[2]);
		if (dose_loaded) {
			this->Draw2DFrom3DDouble(m_spFixedDose, m_spMovingDose, PLANE_AXIAL, curPhysPos[0], m_DoseImgFixed[0], m_DoseImgMoving[0]);
			this->Draw2DFrom3DDouble(m_spFixedDose, m_spMovingDose, PLANE_FRONTAL, curPhysPos[1], m_DoseImgFixed[1], m_DoseImgMoving[1]);
			this->Draw2DFrom3DDouble(m_spFixedDose, m_spMovingDose, PLANE_SAGITTAL, curPhysPos[2], m_DoseImgFixed[2], m_DoseImgMoving[2]);
		}
	}
	else
	{
		this->Draw2DFrom3DDouble(m_spFixed, m_spFixed, PLANE_AXIAL, curPhysPos[0], m_YKImgFixed[0], m_YKImgMoving[0]);
		this->Draw2DFrom3DDouble(m_spFixed, m_spFixed, PLANE_FRONTAL, curPhysPos[1], m_YKImgFixed[1], m_YKImgMoving[1]);
		this->Draw2DFrom3DDouble(m_spFixed, m_spFixed, PLANE_SAGITTAL, curPhysPos[2], m_YKImgFixed[2], m_YKImgMoving[2]);
		if (dose_loaded) {
			this->Draw2DFrom3DDouble(m_spFixedDose, m_spFixedDose, PLANE_AXIAL, curPhysPos[0], m_DoseImgFixed[0], m_DoseImgMoving[0]);
			this->Draw2DFrom3DDouble(m_spFixedDose, m_spFixedDose, PLANE_FRONTAL, curPhysPos[1], m_DoseImgFixed[1], m_DoseImgMoving[1]);
			this->Draw2DFrom3DDouble(m_spFixedDose, m_spFixedDose, PLANE_SAGITTAL, curPhysPos[2], m_DoseImgFixed[2], m_DoseImgMoving[2]);
		}
	}

	//Update position lineEdit
	QString strPos1, strPos2, strPos3;
	strPos1.sprintf("%3.1f", curPhysPos[0]);
	strPos2.sprintf("%3.1f", curPhysPos[1]);
	strPos3.sprintf("%3.1f", curPhysPos[2]);

	ui.lineEditCurPosX->setText(strPos3);
	ui.lineEditCurPosY->setText(strPos2);
	ui.lineEditCurPosZ->setText(strPos1);
	////Update Origin text box
	UShortImageType::PointType imgOriginFixed = m_spFixed->GetOrigin();
	QString strOriFixed;
	strOriFixed.sprintf("%3.4f, %3.4f, %3.4f", imgOriginFixed[0], imgOriginFixed[1], imgOriginFixed[2]);
	ui.lineEditOriginFixed->setText(strOriFixed);

	if (m_spMoving)
	{
		UShortImageType::PointType imgOriginMoving = m_spMoving->GetOrigin();
		QString strOriMoving;
		strOriMoving.sprintf("%3.4f, %3.4f, %3.4f", imgOriginMoving[0], imgOriginMoving[1], imgOriginMoving[2]);
		ui.lineEditOriginMoving->setText(strOriMoving);
	}
	SLT_DrawImageInFixedSlice();
}
//Display is not included here
void gPMC_UI::whenFixedImgLoaded()
{
	if (!m_spFixed)
		return;

	UShortImageType::SizeType imgSize = m_spFixed->GetRequestedRegion().GetSize(); //1016x1016 x z
	UShortImageType::PointType imgOrigin = m_spFixed->GetOrigin();
	UShortImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();

	//to avoid first unnecessary action.
	disconnect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	disconnect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	disconnect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));

	disconnect(ui.sliderFixedW, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));
	disconnect(ui.sliderFixedL, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));
	disconnect(ui.sliderMovingW, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));
	disconnect(ui.sliderMovingL, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));

	initOverlapWndSize();

	ui.sliderPosDisp1->setMinimum(0);
	ui.sliderPosDisp1->setMaximum(imgSize[2] - 1);
	int curPosZ = imgSize[2] / 2;
	ui.sliderPosDisp1->setValue(curPosZ);

	ui.sliderPosDisp2->setMinimum(0);
	ui.sliderPosDisp2->setMaximum(imgSize[1] - 1);
	int curPosY = imgSize[1] / 2;
	ui.sliderPosDisp2->setValue(curPosY);

	ui.sliderPosDisp3->setMinimum(0);
	ui.sliderPosDisp3->setMaximum(imgSize[0] - 1);
	int curPosX = imgSize[0] / 2;
	ui.sliderPosDisp3->setValue(curPosX);

	m_YKDisp[0].SetSplitCenter(QPoint(imgSize[0] / 2, imgSize[1] / 2));
	m_YKDisp[1].SetSplitCenter(QPoint(imgSize[0] / 2, imgSize[2] / 2));
	m_YKDisp[2].SetSplitCenter(QPoint(imgSize[1] / 2, imgSize[2] / 2));

	int iSliderW = 2000;
	int iSliderL = 1024;

	ui.sliderFixedW->setValue(iSliderW);
	ui.sliderMovingW->setValue(iSliderW);

	ui.sliderFixedL->setValue(iSliderL);
	ui.sliderMovingL->setValue(iSliderL);

	connect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	connect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	connect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));

	connect(ui.sliderFixedW, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));
	connect(ui.sliderFixedL, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));
	connect(ui.sliderMovingW, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));
	connect(ui.sliderMovingL, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageInFixedSlice()));
}

void gPMC_UI::SLT_DrawImageInFixedSlice()//Display Swap here!
{
	//Constitute m_YKDisp from Fixed and Moving

	if (ui.checkBoxDrawSplit->isChecked())
	{
		for (int i = 0; i < 3; i++)
		{
			int idxAdd = m_enViewArrange;//m_iViewArrange = 0,1,2
			if (idxAdd + i >= 3)
				idxAdd = idxAdd - 3;

			m_YKDisp[i].SetSpacing(m_YKImgFixed[i + idxAdd].m_fSpacingX, m_YKImgFixed[i + idxAdd].m_fSpacingY);

			m_YKDisp[i].SetSplitOption(PRI_LEFT_TOP);
			if (!m_YKDisp[i].ConstituteFromTwo(m_YKImgFixed[i + idxAdd], m_YKImgMoving[i + idxAdd]))
				std::cout << "Image error " << i + 1 << " th view" << endl;
		}
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			int addedViewIdx = m_enViewArrange;
			if (i + addedViewIdx >= 3)
				addedViewIdx = addedViewIdx - 3;

			m_YKDisp[i].CloneImage(m_YKImgFixed[i + addedViewIdx]);
		}
	}

	// For dose overlay
	if (dose_loaded) {
		if (ui.checkBoxDrawSplit->isChecked())
		{
			for (int i = 0; i < 3; i++)
			{
				int idxAdd = m_enViewArrange;//m_iViewArrange = 0,1,2
				if (idxAdd + i >= 3)
					idxAdd = idxAdd - 3;

				m_AGDisp_Overlay[i].SetSpacing(m_DoseImgFixed[i + idxAdd].m_fSpacingX, m_DoseImgFixed[i + idxAdd].m_fSpacingY);

				m_AGDisp_Overlay[i].SetSplitOption(PRI_LEFT_TOP);
				if (!m_AGDisp_Overlay[i].ConstituteFromTwo(m_DoseImgFixed[i + idxAdd], m_DoseImgMoving[i + idxAdd]))
					std::cout << "Dose Image error " << i + 1 << " th view" << endl;
			}
		}
		else
		{
			for (int i = 0; i < 3; i++)
			{
				int addedViewIdx = m_enViewArrange;
				if (i + addedViewIdx >= 3)
					addedViewIdx = addedViewIdx - 3;

				m_AGDisp_Overlay[i].CloneImage(m_DoseImgFixed[i + addedViewIdx]);
			}
		}
	}

	int sliderW1 = ui.sliderFixedW->value();
	int sliderW2 = ui.sliderMovingW->value();

	int sliderL1 = ui.sliderFixedL->value();
	int sliderL2 = ui.sliderMovingL->value();

	m_YKDisp[0].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
	m_YKDisp[1].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
	m_YKDisp[2].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);

	ui.labelOverlapWnd1->SetBaseImage(&m_YKDisp[0]);
	ui.labelOverlapWnd2->SetBaseImage(&m_YKDisp[1]);
	ui.labelOverlapWnd3->SetBaseImage(&m_YKDisp[2]);

	// here gPMC results could be checked for and displayed, possibly with modification to the qyklabel class /AGA 02/08/2017
	if (dose_loaded) {
		m_AGDisp_Overlay[0].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
		m_AGDisp_Overlay[1].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
		m_AGDisp_Overlay[2].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);

		ui.labelOverlapWnd1->SetOverlayImage(&m_AGDisp_Overlay[0]);
		ui.labelOverlapWnd2->SetOverlayImage(&m_AGDisp_Overlay[1]);
		ui.labelOverlapWnd3->SetOverlayImage(&m_AGDisp_Overlay[2]);
	}

	ui.labelOverlapWnd1->update();
	ui.labelOverlapWnd2->update();
	ui.labelOverlapWnd3->update();
}

void gPMC_UI::SLT_UpdateSplit1()//Mouse Move event
{
	int idx = 0;
	UpdateSplit(idx, ui.labelOverlapWnd1);
}
void gPMC_UI::SLT_UpdateSplit2()
{
	int idx = 1;
	UpdateSplit(idx, ui.labelOverlapWnd2);
}
void gPMC_UI::SLT_UpdateSplit3()
{
	int idx = 2;
	UpdateSplit(idx, ui.labelOverlapWnd3);
}

void gPMC_UI::UpdateSplit(int viewIdx, qyklabel* pOverlapWnd)
{
	int idx = viewIdx;

	if (pOverlapWnd == NULL)
		return;

	if (m_YKDisp[idx].IsEmpty())
		return;

	if (!m_bPressedLeft[idx] && !m_bPressedRight[idx])
		return;

	double dspWidth = pOverlapWnd->width();
	double dspHeight = pOverlapWnd->height();

	int dataWidth = m_YKDisp[idx].m_iWidth;
	int dataHeight = m_YKDisp[idx].m_iHeight;
	if (dataWidth*dataHeight == 0)
		return;

	int dataX = pOverlapWnd->GetDataPtFromMousePos().x();
	int dataY = pOverlapWnd->GetDataPtFromMousePos().y();

	//only works when while the left Mouse is being clicked
	if (m_bPressedLeft[idx])
	{
		m_YKDisp[idx].SetSplitCenter(QPoint(dataX, dataY));
		SLT_DrawImageInFixedSlice();
	}
	else if (m_bPressedRight[idx] && ui.checkBoxPan->isChecked())
	{
		////Update offset information of dispImage

		//GetOriginalDataPos (PanStart)
		//offset should be 0.. only relative distance matters. offset is in realtime changing
		QPoint ptDataPanStartRel = pOverlapWnd->View2DataExt(m_ptPanStart, dspWidth,
			dspHeight, dataWidth, dataHeight, QPoint(0, 0), m_YKDisp[idx].m_fZoom);

		QPoint ptDataPanEndRel = pOverlapWnd->View2DataExt(QPoint(pOverlapWnd->x, pOverlapWnd->y), dspWidth,
			dspHeight, dataWidth, dataHeight, QPoint(0, 0), m_YKDisp[idx].m_fZoom);

		int curOffsetX = ptDataPanEndRel.x() - ptDataPanStartRel.x();
		int curOffsetY = ptDataPanEndRel.y() - ptDataPanStartRel.y();

		int prevOffsetX = m_ptTmpOriginalDataOffset.x();
		int prevOffsetY = m_ptTmpOriginalDataOffset.y();

		m_YKDisp[idx].SetOffset(prevOffsetX - curOffsetX, prevOffsetY - curOffsetY);

		SLT_DrawImageInFixedSlice();
	}
	else if (m_bPressedRight[idx] && !ui.checkBoxPan->isChecked()) //Window Level
	{
		double wWidth = 2.0;
		double wLevel = 2.0;

		int iAddedWidth = (int)((m_ptWindowLevelStart.y() - pOverlapWnd->y)*wWidth);
		int iAddedLevel = (int)((pOverlapWnd->x - m_ptWindowLevelStart.x())*wLevel);

		//Which image is clicked first??
		QPoint crntDataPt = pOverlapWnd->GetDataPtFromViewPt(m_ptWindowLevelStart.x(), m_ptWindowLevelStart.y());

		if (pOverlapWnd->m_pYK16Image != NULL)
		{
			if (m_YKDisp[idx].isPtInFirstImage(crntDataPt.x(), crntDataPt.y()))
			{
				//ui.sliderFixedW->setValue(ui.sliderFixedW->value() + iAddedWidth); //SLT_DrawImageInFixedSlice will be called
				//ui.sliderFixedL->setValue(ui.sliderFixedL->value() + iAddedLevel);
				ui.sliderFixedW->setValue(m_iTmpOriginalW + iAddedWidth); //SLT_DrawImageInFixedSlice will be called
				ui.sliderFixedL->setValue(m_iTmpOriginalL + iAddedLevel);
			}
			else
			{
				ui.sliderMovingW->setValue(m_iTmpOriginalW + iAddedWidth);
				ui.sliderMovingL->setValue(m_iTmpOriginalL + iAddedLevel);
			}
		}
	}
}

//Slide change by scrolling
void gPMC_UI::SLT_MouseWheelUpdate1()
{
	if (ui.checkBoxZoom->isChecked())
	{
		double oldZoom = ui.labelOverlapWnd1->m_pYK16Image->m_fZoom;

		double fWeighting = 0.2;

		ui.labelOverlapWnd1->m_pYK16Image->SetZoom(oldZoom + ui.labelOverlapWnd1->m_iMouseWheelDelta * fWeighting);
		this->SLT_DrawImageInFixedSlice();
	}
	else
		ui.sliderPosDisp1->setValue(ui.sliderPosDisp1->value() + ui.labelOverlapWnd1->m_iMouseWheelDelta);
}
void gPMC_UI::SLT_MouseWheelUpdate2()
{
	if (ui.checkBoxZoom->isChecked())
	{
		double oldZoom = ui.labelOverlapWnd2->m_pYK16Image->m_fZoom;
		double fWeighting = 0.2;

		ui.labelOverlapWnd2->m_pYK16Image->SetZoom(oldZoom + ui.labelOverlapWnd2->m_iMouseWheelDelta * fWeighting);
		this->SLT_DrawImageInFixedSlice();
	}
	else
		ui.sliderPosDisp2->setValue(ui.sliderPosDisp2->value() + ui.labelOverlapWnd2->m_iMouseWheelDelta);
}
void gPMC_UI::SLT_MouseWheelUpdate3()
{
	if (ui.checkBoxZoom->isChecked())
	{
		double oldZoom = ui.labelOverlapWnd3->m_pYK16Image->m_fZoom;
		double fWeighting = 0.2;

		ui.labelOverlapWnd3->m_pYK16Image->SetZoom(oldZoom + ui.labelOverlapWnd3->m_iMouseWheelDelta * fWeighting);
		this->SLT_DrawImageInFixedSlice();
	}
	else
		ui.sliderPosDisp3->setValue(ui.sliderPosDisp3->value() + ui.labelOverlapWnd3->m_iMouseWheelDelta);
}

//release everything
void gPMC_UI::SLT_CancelMouseAction()
{
	for (int i = 0; i < 3; i++)
	{
		m_bPressedLeft[i] = false;//Left Mouse Pressed but not released
		m_bPressedRight[i] = false;
	}
}

void gPMC_UI::SLT_MousePressedLeft1()
{
	m_bPressedLeft[0] = true;
}
void gPMC_UI::SLT_MousePressedLeft2()
{
	m_bPressedLeft[1] = true;
}
void gPMC_UI::SLT_MousePressedLeft3()
{
	m_bPressedLeft[2] = true;
}

void gPMC_UI::MousePressedRight(int wndIdx, qyklabel* pWnd)
{
	m_bPressedRight[wndIdx] = true;

	if (ui.checkBoxPan->isChecked())
	{
		m_ptPanStart.setX(pWnd->x);
		m_ptPanStart.setY(pWnd->y);

		m_ptTmpOriginalDataOffset.setX(m_YKDisp[wndIdx].m_iOffsetX);
		m_ptTmpOriginalDataOffset.setY(m_YKDisp[wndIdx].m_iOffsetY);
	}
	else
	{
		m_ptWindowLevelStart.setX(pWnd->x);
		m_ptWindowLevelStart.setY(pWnd->y);

		QPoint crntDataPt = pWnd->GetDataPtFromViewPt(m_ptWindowLevelStart.x(), m_ptWindowLevelStart.y());

		if (m_YKDisp[wndIdx].isPtInFirstImage(crntDataPt.x(), crntDataPt.y()))
		{
			m_iTmpOriginalL = ui.sliderFixedL->value();
			m_iTmpOriginalW = ui.sliderFixedW->value();
		}
		else
		{
			m_iTmpOriginalL = ui.sliderMovingL->value();
			m_iTmpOriginalW = ui.sliderMovingW->value();
		}
	}
}
void gPMC_UI::SLT_MousePressedRight1()
{
	int idx = 0;
	MousePressedRight(idx, ui.labelOverlapWnd1);
}
void gPMC_UI::SLT_MousePressedRight2()
{
	int idx = 1;
	MousePressedRight(idx, ui.labelOverlapWnd2);
}
void gPMC_UI::SLT_MousePressedRight3()
{
	int idx = 2;
	MousePressedRight(idx, ui.labelOverlapWnd3);
}

void gPMC_UI::SLT_MouseReleasedLeft1()
{
	m_bPressedLeft[0] = false;
}
void gPMC_UI::SLT_MouseReleasedLeft2()
{
	m_bPressedLeft[1] = false;
}
void gPMC_UI::SLT_MouseReleasedLeft3()
{
	m_bPressedLeft[2] = false;
}

void gPMC_UI::SLT_MouseReleasedRight1()
{
	m_bPressedRight[0] = false;
}
void gPMC_UI::SLT_MouseReleasedRight2()
{
	m_bPressedRight[1] = false;
}
void gPMC_UI::SLT_MouseReleasedRight3()
{
	m_bPressedRight[2] = false;
}

void gPMC_UI::SLT_ChangeView()
{
	m_enViewArrange = (m_enViewArrange + 1) % 3;

	YK16GrayImage tmpBufYK[3];

	for (int i = 0; i < 3; i++)
	{
		tmpBufYK[i].CloneImage(m_YKDisp[i]);
	}

	for (int i = 0; i < 3; i++)
	{
		int nextIdx = (i + 1) % 3;
		m_YKDisp[i].CloneImage(tmpBufYK[nextIdx]);
	}
	// /*if (m_enViewArrange > 2)
	//m_enViewArrange = m_enViewArrange-3;*/

	// //Split center should be transfered as well.
	// for (int i = 0 ; i<3 ; i++)
	// {
	//int nextIdx = (i+1)%3;
	///*if (nextIdx >= 3)
	//  nextIdx = nextIdx -3;*/
	//m_YKDisp[i].SetSplitCenter(m_YKDisp[i+1].m_ptSplitCenter);
	////no need of data copy, but other display information (zoom and others should be copied here.
	// }
	initOverlapWndSize();
	shiftSliceSlider();
	updateSliceLabel();

	SLT_DrawImageWhenSliceChange(); //only image data will be updated. Zoom and other things are not.
	//SLT_DrawImageInFixedSlice();
}

void gPMC_UI::initOverlapWndSize()
{
	ui.labelOverlapWnd1->setFixedWidth(DEFAULT_LABEL_SIZE1);
	ui.labelOverlapWnd1->setFixedHeight(DEFAULT_LABEL_SIZE1);

	ui.labelOverlapWnd2->setFixedWidth(DEFAULT_LABEL_SIZE2);
	ui.labelOverlapWnd2->setFixedHeight(DEFAULT_LABEL_SIZE2);

	ui.labelOverlapWnd3->setFixedWidth(DEFAULT_LABEL_SIZE2);
	ui.labelOverlapWnd3->setFixedHeight(DEFAULT_LABEL_SIZE2);
}

void gPMC_UI::shiftSliceSlider() //shift one slice slider information
{
	//to avoid first unnecessary action.
	disconnect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	disconnect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	disconnect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));

	int crntMin[3];
	int crntMax[3];
	int crntValue[3];

	int newMin[3];
	int newMax[3];
	int newValue[3];

	crntMin[0] = ui.sliderPosDisp1->minimum();
	crntMin[1] = ui.sliderPosDisp2->minimum();
	crntMin[2] = ui.sliderPosDisp3->minimum();

	crntMax[0] = ui.sliderPosDisp1->maximum();
	crntMax[1] = ui.sliderPosDisp2->maximum();
	crntMax[2] = ui.sliderPosDisp3->maximum();

	crntValue[0] = ui.sliderPosDisp1->value();
	crntValue[1] = ui.sliderPosDisp2->value();
	crntValue[2] = ui.sliderPosDisp3->value();

	for (int i = 0; i < 3; i++)
	{
		int newIdx = (i + 1) % 3;
		newMin[i] = crntMin[newIdx];
		newMax[i] = crntMax[newIdx];
		newValue[i] = crntValue[newIdx];
	}

	ui.sliderPosDisp1->setMinimum(newMin[0]);
	ui.sliderPosDisp2->setMinimum(newMin[1]);
	ui.sliderPosDisp3->setMinimum(newMin[2]);

	ui.sliderPosDisp1->setMaximum(newMax[0]);
	ui.sliderPosDisp2->setMaximum(newMax[1]);
	ui.sliderPosDisp3->setMaximum(newMax[2]);

	ui.sliderPosDisp1->setValue(newValue[0]);
	ui.sliderPosDisp2->setValue(newValue[1]);
	ui.sliderPosDisp3->setValue(newValue[2]);

	connect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	connect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
	connect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this, SLOT(SLT_DrawImageWhenSliceChange()));
}

void gPMC_UI::updateSliceLabel()
{
	switch (m_enViewArrange)
	{
	case AXIAL_FRONTAL_SAGITTAL:
		ui.labelDisp1->setText("AXIAL");
		ui.labelDisp2->setText("FRONTAL");
		ui.labelDisp3->setText("SAGITTAL");
		break;

	case FRONTAL_SAGITTAL_AXIAL:
		ui.labelDisp1->setText("FRONTAL");
		ui.labelDisp2->setText("SAGITTAL");
		ui.labelDisp3->setText("AXIAL");
		break;

	case SAGITTAL_AXIAL_FRONTAL:
		ui.labelDisp1->setText("SAGITTAL");
		ui.labelDisp2->setText("AXIAL");
		ui.labelDisp3->setText("FRONTAL");
		break;
	}
}

void gPMC_UI::SLT_RestoreImageSingle()
{
	int mainWndIdx = 0;
	m_YKDisp[mainWndIdx].SetZoom(1.0);
	m_YKDisp[mainWndIdx].SetOffset(0, 0);

	SLT_DrawImageInFixedSlice();
}

void gPMC_UI::SLT_RestoreImageAll()
{
	for (int i = 0; i < 3; i++)
	{
		m_YKDisp[i].SetZoom(1.0);
		m_YKDisp[i].SetOffset(0, 0);

		SLT_DrawImageInFixedSlice();
	}
}

void gPMC_UI::Draw2DFrom3DDouble(UShortImageType::Pointer& spFixedImg, UShortImageType::Pointer& spMovingImg, enPLANE enPlane, double pos, YK16GrayImage& YKFixed, YK16GrayImage& YKMoving)
{
	if (!spFixedImg || !spMovingImg)
		return;

	itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(spFixedImg, spFixedImg->GetRequestedRegion());

	UShortImageType::SizeType imgSize = spFixedImg->GetRequestedRegion().GetSize(); //1016x1016 x z
	UShortImageType::PointType imgOrigin = spFixedImg->GetOrigin();
	UShortImageType::SpacingType imgSpacing = spFixedImg->GetSpacing();

	int width = 0;
	int height = 0;
	int iReqSlice = 0;
	int iCntSlice = 0;

	//For moving image
	typedef itk::ResampleImageFilter<UShortImageType, UShortImageType> ResampleFilterType;
	ResampleFilterType::Pointer filter = ResampleFilterType::New();

	filter->SetInput(spMovingImg);

	typedef itk::AffineTransform< double, 3 > TransformType;
	TransformType::Pointer transform = TransformType::New();
	filter->SetTransform(transform);

	typedef itk::NearestNeighborInterpolateImageFunction<UShortImageType, double > InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	filter->SetInterpolator(interpolator);
	filter->SetDefaultPixelValue(0);

	UShortImageType::DirectionType direction;
	direction.SetIdentity();
	filter->SetOutputDirection(direction);

	//ResampledImgType2D::SizeType outSize;

	UShortImageType::SpacingType movingSpacing = imgSpacing;
	UShortImageType::PointType movingOrigin = imgOrigin;
	UShortImageType::SizeType movingSize = imgSize;

	switch (enPlane)
	{
	case PLANE_AXIAL:
		width = imgSize[0];
		height = imgSize[1];
		iCntSlice = imgSize[2];
		iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
		it.SetFirstDirection(0); //x?
		it.SetSecondDirection(1); //y?

		movingSpacing[2] = 1.0;
		movingOrigin[2] = pos;
		movingSize[2] = 1;
		//Resample Here! make corresponding 2D image for Moving image
		YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);

		break;
	case PLANE_FRONTAL:
		width = imgSize[0];
		height = imgSize[2];
		iCntSlice = imgSize[1];
		iReqSlice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
		it.SetFirstDirection(0); //x?
		it.SetSecondDirection(2); //y?

		movingSpacing[1] = 1.0;
		movingOrigin[1] = pos;
		movingSize[1] = 1;

		YKFixed.SetSpacing(imgSpacing[0], imgSpacing[2]);
		break;
	case PLANE_SAGITTAL:
		width = imgSize[1];
		height = imgSize[2];
		iCntSlice = imgSize[0];
		iReqSlice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
		it.SetFirstDirection(1); //x?
		it.SetSecondDirection(2); //y?

		movingSpacing[0] = 1.0;
		movingOrigin[0] = pos;
		movingSize[0] = 1;

		YKFixed.SetSpacing(imgSpacing[1], imgSpacing[2]);
		break;

	default:
		cout << "default should not passed by" << endl;
		width = imgSize[0];
		height = imgSize[1];
		iCntSlice = imgSize[2];
		iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
		it.SetFirstDirection(0); //x?
		it.SetSecondDirection(1); //y?
		YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);
		break;
	}

	filter->SetOutputSpacing(movingSpacing);
	filter->SetOutputOrigin(movingOrigin);
	filter->SetSize(movingSize);
	filter->Update();

	YKFixed.CreateImage(width, height, 0);
	//cout << "Before MovingImg Creation " << endl;

	YKMoving.CreateImage(width, height, 0);//exactly same dimension

	itk::ImageRegionConstIterator<UShortImageType> itMoving(filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
	int cnt = 0;

	//this simple code will cause flip of the image in frontal and sagittal image
	for (itMoving.GoToBegin(); !itMoving.IsAtEnd(); ++itMoving)
	{
		YKMoving.m_pData[cnt] = itMoving.Get();
		cnt++;
	}
	if (enPlane != PLANE_AXIAL)
		YKMoving.EditImage_Flip();

	it.GoToBegin();

	int iNumSlice = 0;
	int iNumWidth = 0;
	int iNumHeight = 0;

	if (iReqSlice < 0 || iReqSlice >= iCntSlice)
		return;

	while (!it.IsAtEnd())
	{
		if (iNumSlice == iReqSlice)
		{
			iNumHeight = 0;

			while (!it.IsAtEndOfSlice())
			{
				iNumWidth = 0;
				while (!it.IsAtEndOfLine())
				{
					UShortImageType::PixelType fixedImgVal = it.Get();
					UShortImageType::IndexType pixelIdxFixed;
					UShortImageType::PointType pixelPhysPt;
					pixelIdxFixed = it.GetIndex();

					if (enPlane == PLANE_AXIAL)
						YKFixed.m_pData[iNumWidth + width*iNumHeight] = fixedImgVal;
					else
						YKFixed.m_pData[iNumWidth + width*(height - iNumHeight - 1)] = fixedImgVal;

					++it;
					iNumWidth++;
				}
				it.NextLine();
				iNumHeight++;
			}
			break;
		}
		it.NextSlice();
		iNumSlice++;
	}
}

// Actually just an overload of the above function
void gPMC_UI::Draw2DFrom3DDouble(UShortImageType::Pointer& spFixedImg, UShortImageType::Pointer& spMovingImg, enPLANE enPlane, double pos, AG17RGBAImage& YKFixed, AG17RGBAImage& YKMoving)
{
	if (!spFixedImg || !spMovingImg)
		return;

	itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(spFixedImg, spFixedImg->GetRequestedRegion());

	UShortImageType::SizeType imgSize = spFixedImg->GetRequestedRegion().GetSize(); //1016x1016 x z
	UShortImageType::PointType imgOrigin = spFixedImg->GetOrigin();
	UShortImageType::SpacingType imgSpacing = spFixedImg->GetSpacing();

	int width = 0;
	int height = 0;
	int iReqSlice = 0;
	int iCntSlice = 0;

	//For moving image
	typedef itk::ResampleImageFilter<UShortImageType, UShortImageType> ResampleFilterType;
	ResampleFilterType::Pointer filter = ResampleFilterType::New();

	filter->SetInput(spMovingImg);

	typedef itk::AffineTransform< double, 3 > TransformType;
	TransformType::Pointer transform = TransformType::New();
	filter->SetTransform(transform);

	typedef itk::NearestNeighborInterpolateImageFunction<UShortImageType, double > InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	filter->SetInterpolator(interpolator);
	filter->SetDefaultPixelValue(0);

	UShortImageType::DirectionType direction;
	direction.SetIdentity();
	filter->SetOutputDirection(direction);

	//ResampledImgType2D::SizeType outSize;

	UShortImageType::SpacingType movingSpacing = imgSpacing;
	UShortImageType::PointType movingOrigin = imgOrigin;
	UShortImageType::SizeType movingSize = imgSize;

	switch (enPlane)
	{
	case PLANE_AXIAL:
		width = imgSize[0];
		height = imgSize[1];
		iCntSlice = imgSize[2];
		iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
		it.SetFirstDirection(0); //x?
		it.SetSecondDirection(1); //y?

		movingSpacing[2] = 1.0;
		movingOrigin[2] = pos;
		movingSize[2] = 1;
		//Resample Here! make corresponding 2D image for Moving image
		YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);

		break;
	case PLANE_FRONTAL:
		width = imgSize[0];
		height = imgSize[2];
		iCntSlice = imgSize[1];
		iReqSlice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
		it.SetFirstDirection(0); //x?
		it.SetSecondDirection(2); //y?

		movingSpacing[1] = 1.0;
		movingOrigin[1] = pos;
		movingSize[1] = 1;

		YKFixed.SetSpacing(imgSpacing[0], imgSpacing[2]);
		break;
	case PLANE_SAGITTAL:
		width = imgSize[1];
		height = imgSize[2];
		iCntSlice = imgSize[0];
		iReqSlice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
		it.SetFirstDirection(1); //x?
		it.SetSecondDirection(2); //y?

		movingSpacing[0] = 1.0;
		movingOrigin[0] = pos;
		movingSize[0] = 1;

		YKFixed.SetSpacing(imgSpacing[1], imgSpacing[2]);
		break;

	default:
		cout << "default should not passed by" << endl;
		width = imgSize[0];
		height = imgSize[1];
		iCntSlice = imgSize[2];
		iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
		it.SetFirstDirection(0); //x?
		it.SetSecondDirection(1); //y?
		YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);
		break;
	}

	filter->SetOutputSpacing(movingSpacing);
	filter->SetOutputOrigin(movingOrigin);
	filter->SetSize(movingSize);
	filter->Update();

	YKFixed.CreateImage(width, height, 0);
	YKMoving.CreateImage(width, height, 0);//exactly same dimension

	itk::ImageRegionConstIterator<UShortImageType> itMoving(filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
	int cnt = 0;

	//this simple code will cause flip of the image in frontal and sagittal image
	for (itMoving.GoToBegin(); !itMoving.IsAtEnd(); ++itMoving)
	{
		YKMoving.m_pData[cnt] = itMoving.Get();
		cnt++;
	}
	if (enPlane != PLANE_AXIAL)
	{
		YKMoving.EditImage_Flip();
	}

	it.GoToBegin();

	int iNumSlice = 0;
	int iNumWidth = 0;
	int iNumHeight = 0;

	if (iReqSlice < 0 || iReqSlice >= iCntSlice)
		return;

	while (!it.IsAtEnd())
	{
		if (iNumSlice == iReqSlice)
		{
			iNumHeight = 0;

			while (!it.IsAtEndOfSlice())
			{
				iNumWidth = 0;
				while (!it.IsAtEndOfLine())
				{
					UShortImageType::PixelType fixedImgVal = it.Get();

					if (enPlane == PLANE_AXIAL)
						YKFixed.m_pData[iNumWidth + width*iNumHeight] = fixedImgVal;
					else
						YKFixed.m_pData[iNumWidth + width*(height - iNumHeight - 1)] = fixedImgVal;

					++it;
					iNumWidth++;
				}
				it.NextLine();
				iNumHeight++;
			}
			break;
		}
		it.NextSlice();
		iNumSlice++;
	}
}

bool gPMC_UI::eventFilter(QObject *target, QEvent *event)
{
	if (event->type() == QEvent::Paint) //eventy->type() = 12 if paint event
		return false;

	return QWidget::eventFilter(target, event);//This is mandatory to deliver the event signal to the children components
}

void gPMC_UI::childEvent(QChildEvent *event)
{
	//cout << "childEvent Called" << endl;
	if (event->added()) {
		event->child()->installEventFilter(this);
	}
}

QString get_output_options(UShortImageType::Pointer m_spFixed){
	QString str_fixedOrigin = QString("%1,%2,%3") // done per image because CT might be different from reconstructed CBCT
		.arg(m_spFixed->GetOrigin()[0])
		.arg(m_spFixed->GetOrigin()[1])
		.arg(m_spFixed->GetOrigin()[2]);
	QString str_fixedDimension = QString("%1,%2,%3")
		.arg(m_spFixed->GetBufferedRegion().GetSize()[0])
		.arg(m_spFixed->GetBufferedRegion().GetSize()[1])
		.arg(m_spFixed->GetBufferedRegion().GetSize()[2]);
	QString str_fixedSpacing = QString("%1,%2,%3")
		.arg(m_spFixed->GetSpacing()[0])
		.arg(m_spFixed->GetSpacing()[1])
		.arg(m_spFixed->GetSpacing()[2]);
	QString str_fixedDirection = QString("%1,%2,%3,%4,%5,%6,%7,%8,%9")
		.arg(m_spFixed->GetDirection()[0][0])
		.arg(m_spFixed->GetDirection()[0][1])
		.arg(m_spFixed->GetDirection()[0][2])
		.arg(m_spFixed->GetDirection()[1][0])
		.arg(m_spFixed->GetDirection()[1][1])
		.arg(m_spFixed->GetDirection()[1][2])
		.arg(m_spFixed->GetDirection()[2][0])
		.arg(m_spFixed->GetDirection()[2][1])
		.arg(m_spFixed->GetDirection()[2][2]);

	return QString(" --origin=%1 --spacing=%2 --dimension=%3 --direction=%4")
		.arg(str_fixedOrigin)
		.arg(str_fixedSpacing)
		.arg(str_fixedDimension)
		.arg(str_fixedDirection);
}

inline void create_and_exec_gPMC_cmd(QString plan_path, QString fixed_common_str, QString moving_common_str){
	QString gPMC_command_str = QString("%1 --plan=\"%2\"")
		.arg(fixed_common_str).arg(plan_path);

	std::cout << gPMC_command_str.toStdString() << std::endl;
	if (QProcess::execute(gPMC_command_str) < 0)
		std::cerr << "Failed to run (fixed mc recalc)" << std::endl;

	if (moving_common_str != "")
	{
		gPMC_command_str = QString("%1 --plan=\"%2\"")
			.arg(moving_common_str).arg(plan_path);

		if (QProcess::execute(gPMC_command_str) < 0)
			std::cerr << "Failed to run (moving mc recalc)" << std::endl;
	}
}

UShortImageType::Pointer Read_Add_Normalize_Combine(QStringList image_paths){
	typedef itk::ImageFileReader<FloatImageType> ImageReaderType;
	typedef itk::MinimumMaximumImageCalculator<FloatImageType> MinMaxFindType;
	typedef itk::MultiplyImageFilter<FloatImageType, FloatImageType, UShortImageType> MultiplyImageFilterType;
	typedef itk::AddImageFilter<FloatImageType, FloatImageType> AddImageFilterType;

	ImageReaderType::Pointer DoseReader = ImageReaderType::New();
	DoseReader->SetFileName(image_paths[0].toStdString());
	DoseReader->UpdateOutputInformation();
	DoseReader->Update();
	FloatImageType::Pointer tmp_Img = DoseReader->GetOutput();

	for (int i = 1; i < image_paths.size(); i++){
		ImageReaderType::Pointer DoseReader2 = ImageReaderType::New();
		DoseReader2->SetFileName(image_paths[i].toStdString());
		DoseReader2->UpdateOutputInformation();
		DoseReader2->Update();
		AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
		addFilter->SetInput1(tmp_Img);
		addFilter->SetInput2(DoseReader2->GetOutput());
		addFilter->Update();
		tmp_Img = addFilter->GetOutput(); // outOfPlace because inplace apperently doesn't propagate buffered region properly
	}

	MinMaxFindType::Pointer MinMaxFilter = MinMaxFindType::New();
	MinMaxFilter->SetImage(tmp_Img);
	MinMaxFilter->ComputeMaximum();
	std::cout << "Max [MeV/g/primary]: " << MinMaxFilter->GetMaximum() << " scaled to 65535" << std::endl;
	//Multiply: Scale to USHORT
	MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
	multiplyImageFilter->SetInput(tmp_Img);
	multiplyImageFilter->SetConstant(65535.0f / MinMaxFilter->GetMaximum());

	multiplyImageFilter->Update();
	return multiplyImageFilter->GetOutput();
}

void gPMC_UI::SLT_recalcDose()
{
	if (!m_spFixed)
		return;
	ui.progressBar->setValue(1);
	QString fixed_dcm_dir = "";
	QString moving_dcm_dir = "";
	// Export fixed and moving as DCM
	fixed_dcm_dir = SaveUSHORTAsSHORT_DICOM_gdcmITK(m_spFixed, QString("tmp_"), QString("Fixed"), m_strPathDirDefault).absolutePath();
	if (m_spFixed != m_spMoving)
		moving_dcm_dir = SaveUSHORTAsSHORT_DICOM_gdcmITK(m_spMoving, QString("tmp_"), QString("Moving"), m_strPathDirDefault).absolutePath();

	std::cout << fixed_dcm_dir.toStdString() << std::endl;
	ui.progressBar->setValue(15);

	// std::vector<QString> plan_filepaths;
	// Load dcm rtplan.
	QStringList plan_filepaths = QFileDialog::getOpenFileNames(this, "Open DCMRT Plan file",
		m_strPathDirDefault, "DCMRT Plan (*.dcm)", 0, 0);
	size_t n_plans = (size_t)plan_filepaths.size();
	/* Below MUST be done externally, through command line - due to incompatibilities with compilers and library versions..
	* Pytrip-style but with somewhat better control.. (hacks to private members might be possible as well: http://stackoverflow.com/questions/424104/can-i-access-private-members-from-outside-the-class-without-using-friends)
	*/

	// Create gPMC command line

	QString gPMC_device;
	gPMC_device = "cpu";

	if (ui.radioButton_UseGPU->isChecked())
		gPMC_device = "gpu";
	else if (ui.radioButton_UseACC->isChecked())
		gPMC_device = "acc";

	//"dose2water","letd","dose2medium","fluence"
	QString gPMC_modality;
	gPMC_modality = "dose2medium";

	if (ui.radioButton_calcLET->isChecked())
		gPMC_modality = "letd";
	else if (ui.radioButton_calcD2W->isChecked())
		gPMC_modality = "dose2water";
	else if (ui.radioButton_calcFluence->isChecked())
		gPMC_modality = "fluence";

	QString image_independent_string = QString(" --metric=\"%1\" --hardware=\"%2\" -b %3 --verbose")
		.arg(gPMC_modality)
		.arg(gPMC_device)
		.arg(ui.spinBox_Nsims->value());
	//.arg(n_plans > 1 ? "" : " --verbose");

	QStringList FixedImagePaths;
	QStringList MovingImagePaths;

	for (size_t i = 0; i < n_plans; i++){
		FixedImagePaths.push_back(QString("%1/dose_fixed_%2.mha").arg(fixed_dcm_dir).arg(i + 1));
		MovingImagePaths.push_back(QString("%1/dose_moving_%2.mha").arg(moving_dcm_dir).arg(i + 1));
	}

	ui.progressBar->setValue(20);

	for (size_t i = 0; i < n_plans; i++){
		QString gPMC_command_str = QString("gPMC.exe --dir=\"%1\" --output=\"%2\"%3%4")
			.arg(fixed_dcm_dir)
			.arg(QFileInfo(FixedImagePaths[i]).absoluteFilePath())
			.arg(get_output_options(m_spFixed))
			.arg(image_independent_string);

		QString gPMC_mov_command_str = "";
		if (moving_dcm_dir != "")
		{
			gPMC_mov_command_str = QString("gPMC.exe --dir=\"%1\" --output=\"%2\"%3%4")
				.arg(moving_dcm_dir)
				.arg(QFileInfo(MovingImagePaths[i]).absoluteFilePath())
				.arg(get_output_options(m_spMoving))
				.arg(image_independent_string);
		}

		create_and_exec_gPMC_cmd(plan_filepaths[i], gPMC_command_str, gPMC_mov_command_str);
		ui.progressBar->setValue((int)(20 + (60 * i) / n_plans));
	}
	// Run gPMC externally ^

	// Translate gPMC output (preferably .mha) to ITK image
	QFileInfo finfoDosePath = QFileInfo(FixedImagePaths[0]);
	if (!finfoDosePath.exists())
		return;

	m_spFixedDose = Read_Add_Normalize_Combine(FixedImagePaths);

	if (!m_spFixedDose)
		std::cout << "Dose failed to load for fixed Image!!" << std::endl;
	else
		std::cout << "Dose loaded for fixed Image" << std::endl;

	UShortImageType::SizeType imgDim = m_spFixedDose->GetBufferedRegion().GetSize();
	UShortImageType::SpacingType spacing = m_spFixedDose->GetSpacing();

	std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1] << "	" << imgDim[2] << std::endl;
	std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1] << "	" << spacing[2] << std::endl;

	ui.progressBar->setValue(90);
	if (moving_dcm_dir != "")
	{
		m_spMovingDose = Read_Add_Normalize_Combine(MovingImagePaths);
		std::cout << "Dose loaded for moving Image" << std::endl;
	}
	else
	{
		if (!m_spFixedDose.IsNull())
			m_spMovingDose = m_spFixedDose;
	}
	// Display dose as colorwash on top of fixed and moving in all three plots
	if (!m_spFixedDose && !m_spMovingDose)
		dose_loaded = false;
	else
		dose_loaded = true;

	ui.progressBar->setValue(99);
	SLT_DrawImageWhenSliceChange();

	ui.progressBar->setValue(0);
}

void gPMC_UI::SLT_FixedImageSelected(QString selText)
{
	LoadImgFromComboBox(0, selText); // here, m_spMoving and Fixed images are determined
}

void gPMC_UI::SLT_MovingImageSelected(QString selText)
{
	LoadImgFromComboBox(1, selText);
}

void gPMC_UI::LoadImgFromComboBox(int idx, QString& strSelectedComboTxt)// -->when fixed image loaded will be called here!
{
	//cout << "LoadImgFromComboBox " << "index " << idx << "text " << strSelectedComboTxt.toLocal8Bit().constData() << endl;

	UShortImageType::Pointer spTmpImg;
	if (strSelectedComboTxt.compare(QString("RAW_CBCT"), Qt::CaseSensitive) == 0)
		spTmpImg = m_spRawReconImg;
	else if (strSelectedComboTxt.compare(QString("REF_CT"), Qt::CaseSensitive) == 0)
		spTmpImg = m_spRefCTImg;

	if (!spTmpImg)
		return;

	if (idx == 0)
	{
		m_spFixed = spTmpImg;

		whenFixedImgLoaded();
	}
	else if (idx == 1)
	{
		m_spMoving = spTmpImg;
	}

	SLT_DrawImageWhenSliceChange();
}

//search  for the  main data, if there  is, add  the predefined name to the combobox
void gPMC_UI::UpdateListOfComboBox(int idx)
{
	QComboBox* crntCombo;

	if (idx == 0)
		crntCombo = ui.comboBoxImgFixed;
	else
		crntCombo = ui.comboBoxImgMoving;

	//remove all the list
	crntCombo->clear();

	if (m_spRawReconImg)
		crntCombo->addItem("RAW_CBCT");

	if (m_spRefCTImg)
		crntCombo->addItem("REF_CT");
}

//Bring Focus
void gPMC_UI::SLT_BringFocusToEnableArrow(bool bChecked)
{
	if (bChecked)
		ui.labelDisp1->setFocus(); //if focus is in label, Key Event will trigger 51 (override)
}

void gPMC_UI::SLT_KeyMoving(bool bChecked)//Key Moving check box
{
	ui.comboBoxImgFixed->setDisabled(bChecked);
	ui.comboBoxImgMoving->setDisabled(bChecked);
}

void gPMC_UI::SLTM_LoadDICOMdir()
{
	ui.progressBar->setValue(1);
	QString dirPath = QFileDialog::getExistingDirectory(this, tr("Open Directory"), m_strPathDirDefault, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

	if (dirPath.length() <= 1)
		return;

	ui.progressBar->setValue(20);

	itk::GDCMImageIO::Pointer gdcmIO = itk::GDCMImageIO::New();
	itk::GDCMSeriesFileNames::Pointer inputNames = itk::GDCMSeriesFileNames::New();
	inputNames->SetInputDirectory(dirPath.toStdString());

	typedef itk::ImageSeriesReader< ShortImageType > ReaderType;
	const ReaderType::FileNamesContainer & filenames = inputNames->GetInputFileNames();
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetImageIO(gdcmIO);
	reader->SetFileNames(filenames);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Exception thrown while reading the series" << std::endl;
		std::cerr << excp << std::endl;
		return;
	}

	ui.progressBar->setValue(90);

	typedef itk::MinimumMaximumImageCalculator <ShortImageType> ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(reader->GetOutput());
	imageCalculatorFilter->Compute();

	double minVal = (double)(imageCalculatorFilter->GetMinimum());
	double maxVal = (double)(imageCalculatorFilter->GetMaximum());

	cout << "Current Min and Max Values are	" << minVal << "	" << maxVal << endl;
	ui.progressBar->setValue(95);

	unsigned short outputMinVal = (unsigned short)(minVal + 1024);
	unsigned short outputMaxVal = (unsigned short)(maxVal + 1024);

	typedef itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType> RescaleFilterType;
	RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
	spRescaleFilter->SetInput(reader->GetOutput());
	spRescaleFilter->SetOutputMinimum(outputMinVal);
	spRescaleFilter->SetOutputMaximum(outputMaxVal);
	spRescaleFilter->Update();

	ui.progressBar->setValue(99);

	//m_spRawReconImg = spRescaleFilter->GetOutput();
	m_spRefCTImg = spRescaleFilter->GetOutput();

	m_spFixed = m_spRefCTImg;
	if (m_spMoving.IsNull())
		m_spMoving = m_spFixed;

	UpdateListOfComboBox(0);
	SLT_DrawImageWhenSliceChange();

	ui.progressBar->setValue(0);
}

ScatterDataModifier* setup_CustomWidget(QWidget* widget){
	using namespace QtDataVisualization;
	{
		Q3DScatter *graph = new Q3DScatter();
		QWidget* container = QWidget::createWindowContainer(graph);
		if (!graph->hasContext()) {
			std::cerr << "Couldn't initialize the OpenGL context." << std::endl;
			return NULL;
		}

		QSize screenSize = graph->screen()->size();
		container->setMinimumSize(QSize(screenSize.width() / 2, screenSize.height() / 1.5));
		container->setMaximumSize(screenSize);
		container->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		container->setFocusPolicy(Qt::StrongFocus);

		QHBoxLayout *hLayout = new QHBoxLayout(widget);
		QVBoxLayout *vLayout = new QVBoxLayout();
		hLayout->addWidget(container, 1);
		hLayout->addLayout(vLayout);

		// widget->setWindowTitle(QStringLiteral("A Cosine Wave"));

		QComboBox *themeList = new QComboBox(widget);
		themeList->addItem(QStringLiteral("Qt"));
		themeList->addItem(QStringLiteral("Primary Colors"));
		themeList->addItem(QStringLiteral("Digia"));
		themeList->addItem(QStringLiteral("Stone Moss"));
		themeList->addItem(QStringLiteral("Army Blue"));
		themeList->addItem(QStringLiteral("Retro"));
		themeList->addItem(QStringLiteral("Ebony"));
		themeList->addItem(QStringLiteral("Isabelle"));
		themeList->setCurrentIndex(0);

		QPushButton *labelButton = new QPushButton(widget);
		labelButton->setText(QStringLiteral("Change label style"));

		QCheckBox *smoothCheckBox = new QCheckBox(widget);
		smoothCheckBox->setText(QStringLiteral("Smooth dots"));
		smoothCheckBox->setChecked(false);

		QComboBox *itemStyleList = new QComboBox(widget);
		itemStyleList->addItem(QStringLiteral("Sphere"), int(QAbstract3DSeries::MeshSphere));
		itemStyleList->addItem(QStringLiteral("Cube"), int(QAbstract3DSeries::MeshCube));
		itemStyleList->addItem(QStringLiteral("Minimal"), int(QAbstract3DSeries::MeshMinimal));
		itemStyleList->addItem(QStringLiteral("Point"), int(QAbstract3DSeries::MeshPoint));
		itemStyleList->setCurrentIndex(3);

		QPushButton *cameraButton = new QPushButton(widget);
		cameraButton->setText(QStringLiteral("Change camera preset"));
#ifdef QT_TEST
		QPushButton *itemCountButton = new QPushButton(widget);
		itemCountButton->setText(QStringLiteral("Toggle item count"));
#endif
		QCheckBox *backgroundCheckBox = new QCheckBox(widget);
		backgroundCheckBox->setText(QStringLiteral("Show background"));
		backgroundCheckBox->setChecked(true);

		QCheckBox *gridCheckBox = new QCheckBox(widget);
		gridCheckBox->setText(QStringLiteral("Show grid"));
		gridCheckBox->setChecked(true);

		QComboBox *shadowQuality = new QComboBox(widget);
		shadowQuality->addItem(QStringLiteral("None"));
		shadowQuality->addItem(QStringLiteral("Low"));
		shadowQuality->addItem(QStringLiteral("Medium"));
		shadowQuality->addItem(QStringLiteral("High"));
		shadowQuality->addItem(QStringLiteral("Low Soft"));
		shadowQuality->addItem(QStringLiteral("Medium Soft"));
		shadowQuality->addItem(QStringLiteral("High Soft"));
		shadowQuality->setCurrentIndex(0);

		vLayout->addWidget(labelButton, 0, Qt::AlignTop);
		vLayout->addWidget(cameraButton, 0, Qt::AlignTop);
#ifdef QT_TEST
		vLayout->addWidget(itemCountButton, 0, Qt::AlignTop);
#endif
		vLayout->addWidget(backgroundCheckBox);
		vLayout->addWidget(gridCheckBox);
		vLayout->addWidget(smoothCheckBox, 0, Qt::AlignTop);
		vLayout->addWidget(new QLabel(QStringLiteral("Change dot style")));
		vLayout->addWidget(itemStyleList);
		vLayout->addWidget(new QLabel(QStringLiteral("Change theme")));
		vLayout->addWidget(themeList);
		vLayout->addWidget(new QLabel(QStringLiteral("Adjust shadow quality")));
		vLayout->addWidget(shadowQuality);

		ScatterDataModifier* m_modifier = new ScatterDataModifier(graph);

		QObject::connect(cameraButton, &QPushButton::clicked, m_modifier,
			&ScatterDataModifier::changePresetCamera);
		QObject::connect(labelButton, &QPushButton::clicked, m_modifier,
			&ScatterDataModifier::changeLabelStyle);
#ifdef QT_TEST
		QObject::connect(itemCountButton, &QPushButton::clicked, m_modifier,
			&ScatterDataModifier::toggleItemCount);
#endif
		QObject::connect(backgroundCheckBox, &QCheckBox::stateChanged, m_modifier,
			&ScatterDataModifier::setBackgroundEnabled);
		QObject::connect(gridCheckBox, &QCheckBox::stateChanged, m_modifier,
			&ScatterDataModifier::setGridEnabled);
		QObject::connect(smoothCheckBox, &QCheckBox::stateChanged, m_modifier,
			&ScatterDataModifier::setSmoothDots);

		QObject::connect(m_modifier, &ScatterDataModifier::backgroundEnabledChanged,
			backgroundCheckBox, &QCheckBox::setChecked);
		QObject::connect(m_modifier, &ScatterDataModifier::gridEnabledChanged,
			gridCheckBox, &QCheckBox::setChecked);
		QObject::connect(itemStyleList, SIGNAL(currentIndexChanged(int)), m_modifier,
			SLOT(changeStyle(int)));

		QObject::connect(themeList, SIGNAL(currentIndexChanged(int)), m_modifier,
			SLOT(changeTheme(int)));

		QObject::connect(shadowQuality, SIGNAL(currentIndexChanged(int)), m_modifier,
			SLOT(changeShadowQuality(int)));

		QObject::connect(m_modifier, &ScatterDataModifier::shadowQualityChanged, shadowQuality,
			&QComboBox::setCurrentIndex);
		QObject::connect(graph, &Q3DScatter::shadowQualityChanged, m_modifier,
			&ScatterDataModifier::shadowQualityUpdatedByVisual);

		return m_modifier;
	}
}

void gPMC_UI::SLT_pltBeamGeo(){
	QWidget *widget = new QWidget; // <- If you want a new window instead of embedded
	m_modifier = setup_CustomWidget(widget); //ui.CustomWidget);
	load_dcm_into_ScatterWidget();
	widget->show();
}

void gPMC_UI::load_dcm_into_ScatterWidget(){
	size_t N_dicom = 0;

	QStringList plan_paths = QFileDialog::getOpenFileNames(this, "Open DCMRT Plan file",
		m_strPathDirDefault, "DCMRT Plan (*.dcm)", 0, 0);

	std::vector<std::string> plan_vector;
	for (int i = 0; i < plan_paths.size(); i++)
		plan_vector.push_back(plan_paths[i].toStdString());

	for (size_t i = 0; i < plan_vector.size(); i++)
		N_dicom += GetNFromDicom<10>(plan_vector[i]);

	m_modifier->setItemCount(N_dicom);

	cl_float* T = new cl_float[N_dicom];      //Energy(MeV?)= nominal control point energy
	cl_float3* pos = new cl_float3[N_dicom];  //Position    = Scan Spot Position transformed to gPMC coordinates
	cl_float3* dir = new cl_float3[N_dicom];  //Direction   = vector in unit sphere defined by the gantry and couch angle
	cl_float* weight = new cl_float[N_dicom]; //Weight      = Scan Spot Meterset Weights of control point

	size_t N_result = initSourceFromDicom<10>(plan_vector, T, pos, dir, weight);

	if (N_result != N_dicom)
		std::cout << "\a" << "ONE OF THE COUNTERS ARE WRONG!!" << " Dicom: " << N_dicom << ", Result: " << N_result << std::endl;

	m_modifier->setPosData(pos);
	m_modifier->setDirectionData(dir);
	m_modifier->setWeightData(weight);
	m_modifier->setEnergyData(T);
	m_modifier->callAddData();
}