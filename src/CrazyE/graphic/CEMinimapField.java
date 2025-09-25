package CrazyE.graphic;

import arc.math.geom.Vec2;

public interface CEMinimapField {
    default int color(int i, int j){
        return 0;
    }

    default Vec2 arrowSize(int i, int j){
        return Vec2.ZERO;
    }

    enum FieldType{
        Heatmap,; //Vectorfield;
    }
}
