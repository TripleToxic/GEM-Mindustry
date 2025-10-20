package CrazyE.world.environment;

import CrazyE.graphic.*;
import arc.Events;
import arc.func.Cons;
import arc.math.Mathf;
import arc.math.geom.*;
import arc.struct.Seq;
import arc.util.*;
import mindustry.game.*;
import mindustry.gen.Unit;
import mindustry.io.SaveFileReader.*;
import java.io.*;

import static java.lang.Math.*;
import static mindustry.Vars.*;
import static mindustry.game.EventType.*;
import static mindustry.io.SaveVersion.addCustomChunk;
import static CrazyE.graphic.CEMinimapRenderer.*;

public class Gravity implements CustomChunk, CEMinimapField{
    private static final double cons = cbrt(2),
            w1 = -(cons/(2d - cons)),
            w2 = (1d/(2d - cons));
    private static final float c = 20f * tilesize * tilesize / 60f, maxSpeed = 0.9999f*c, csquared = c * c, G = 0.1f, ms = 4f, stabilizer = 1.1f, timescale = 1/32f,
    c1 = (float)(0.5d*w2),
    c2 = (float)(0.5d*(w1+w2)),
    c3 = c2,
    c4 = c1,

    d1 = (float)w2,
    d2 = (float)w1,
    d3 = (float)w2,
    drawScale = 16f;

    private float dt, ac, b, b_inv;
    private int sizeX, sizeY;
    private GEMgrid PotentialField, vPotentialField, vBuffer, BufferField, BufferField2, CurrentField, EMGravityField;
    private final float[] fm = new float[9];
    private float[] tridParam1, tridParam2;
    private Vec3[] tridResult;

    String debug = "";

    public Gravity(){
        Events.on(WorldLoadEvent.class, e -> {
            reset();
        });

        Events.run(Trigger.beforeGameUpdate, () -> {
            if(state.rules.editor) return;
            CurrentField.clear();
            initCurrent();
        });

        Events.run(Trigger.afterGameUpdate, () -> {
            if(state.rules.editor) return;
            dt = Time.delta*timescale;//Mathf.clamp(Time.delta, 0, 1.5f)*timescale;
            ac = -0.25f*csquared*dt*dt; //in this tridiagonal matrix algorithm, a == c, so I put them together into one name
            b = 1f - 2f*ac;
            b_inv = 1f/b;

            initCurrent();
            CurrentField.scl(0.5f);

            tridParam1[0] = 2*ac*b_inv;

            for(int i = 1; i < tridParam1.length-1; i++){
                tridParam2[i] = 1f/(b - ac*tridParam1[i-1]);
                tridParam1[i] = ac*tridParam2[i];
            }

            tridParam2[tridParam2.length-1] = 1f/(b - 2*ac*tridParam1[tridParam2.length-1]);

            updateFields();
            updateOther();
        });

        Events.on(SaveWriteEvent.class, e -> {
            Log.info(debug);
        });

        addField("Gravity", this);
        addCustomChunk("Gravity", this);
    }

    public void reset(){
        this.sizeX = GEMgrid.sizeX = world.width();
        this.sizeY = GEMgrid.sizeY = world.height();

        int maxSize = max(sizeX, sizeY);

        this.tridParam1 = new float[maxSize];
        this.tridParam2 = new float[maxSize];
        this.tridResult = new Vec3[maxSize];
        for(int i = 0; i < maxSize; i++){
            tridResult[i] = new Vec3();
        }
        PotentialField = new GEMgrid();
        vPotentialField = new GEMgrid();
        vBuffer = new GEMgrid();
        BufferField = new GEMgrid();
        BufferField2 = new GEMgrid();
        CurrentField = new GEMgrid();
        EMGravityField = new GEMgrid();
    }

    public void EMGravity(){
        float Ex, Ey, B;
        for(int y=1; y<sizeY-1; y++){
            for(int x=1; x<sizeX-1; x++){
                Ex = -0.5f*(PotentialField.get(y, x + 1).x - PotentialField.get(y, x - 1).x);
                Ey = -0.5f*(PotentialField.get(y + 1, x).x - PotentialField.get(y - 1, x).x);

                Ex -= vPotentialField.get(y, x).y;
                Ey -= vPotentialField.get(y, x).z;

                B = 0.5f*(PotentialField.get(y, x + 1).z - PotentialField.get(y, x - 1).z)
                        - 0.5f*(PotentialField.get(y + 1, x).y - PotentialField.get(y - 1, x).y);
                EMGravityField.get(y, x).set(Ex, Ey, B);
            }
        }

        for(int y=1; y<sizeY-1; y++){
            Ex = -vPotentialField.get(y, 0).x/c;
            Ey = -0.5f*(PotentialField.get(y + 1, 0).x - PotentialField.get(y - 1, 0).x);

            Ex -= vPotentialField.get(y, 0).y;
            Ey -= vPotentialField.get(y, 0).z;

            B = vPotentialField.get(y, 0).z/c -
                    - 0.5f*(PotentialField.get(y + 1, 0).y - PotentialField.get(y - 1, 0).y);

            EMGravityField.get(y, 0).set(Ex, Ey, B);

            Ex = vPotentialField.get(y, sizeX-1).x/c;
            Ey = -0.5f*(PotentialField.get(y + 1, sizeX-1).x - PotentialField.get(y - 1, sizeX-1).x);

            Ex -= vPotentialField.get(y, sizeX-1).y;
            Ey -= vPotentialField.get(y, sizeX-1).z;

            B = -vPotentialField.get(y, sizeX-1).z/c -
                    - 0.5f*(PotentialField.get(y + 1, sizeX-1).y - PotentialField.get(y - 1, sizeX-1).y);

            EMGravityField.get(y, sizeX-1).set(Ex, Ey, B);
        }

        for(int x=1; x<sizeX-1; x++){
            Ex = -0.5f*(PotentialField.get(0, x + 1).x - PotentialField.get(0, x - 1).x);
            Ey = -vPotentialField.get(0, x).x/c;

            Ex -= vPotentialField.get(0, x).y;
            Ey -= vPotentialField.get(0, x).z;

            B = 0.5f*(PotentialField.get(0, x + 1).z - PotentialField.get(0, x - 1).z)
                    - vPotentialField.get(0, x).y/c;

            EMGravityField.get(0, x).set(Ex, Ey, B);

            Ex = -0.5f*(PotentialField.get(sizeY-1, x + 1).x - PotentialField.get(sizeY-1, x - 1).x);
            Ey = vPotentialField.get(sizeY-1, x).x/c;

            Ex -= vPotentialField.get(sizeY-1, x).y;
            Ey -= -vPotentialField.get(sizeY-1, x).z;

            B = 0.5f*(PotentialField.get(sizeY-1, x + 1).z - PotentialField.get(sizeY-1, x - 1).z)
                     + vPotentialField.get(sizeY-1, x).z/c;

            EMGravityField.get(sizeY-1, x).set(Ex, Ey, B);
        }

        //Bottom-Left
        Ey = Ex = -vPotentialField.get(0, 0).x/c;

        Ex -= vPotentialField.get(0, 0).y;
        Ey -= vPotentialField.get(0, 0).z;

        B = vPotentialField.get(0,0).z
                - vPotentialField.get(0, 0).y;

        EMGravityField.get(0, 0).set(Ex, Ey, B/c);

        //Bottom-Right
        Ex = vPotentialField.get(0, sizeX-1).x/c;
        Ey = -Ex;

        Ex -= vPotentialField.get(0, sizeX-1).y;
        Ey -= vPotentialField.get(0, sizeX-1).z;

        B = -vPotentialField.get(0,sizeX-1).z
                - vPotentialField.get(0, sizeX-1).y;

        EMGravityField.get(0, sizeX-1).set(Ex, Ey, B/c);

        //Top-Left
        Ex = -vPotentialField.get(sizeY-1, 0).x/c;
        Ey = -Ex;

        Ex -= vPotentialField.get(sizeY-1, 0).y;
        Ey -= vPotentialField.get(sizeY-1, 0).z;

        B = vPotentialField.get(sizeY-1,0).z
                + vPotentialField.get(sizeY-1, 0).y;

        EMGravityField.get(sizeY-1, 0).set(Ex, Ey, B/c);

        //Top-Right
        Ey = Ex = vPotentialField.get(sizeY-1, sizeX-1).x/c;

        Ex -= vPotentialField.get(sizeY-1, sizeX-1).y;
        Ey -= vPotentialField.get(sizeY-1, sizeX-1).z;

        B = -vPotentialField.get(sizeY-1, sizeX-1).z
                + vPotentialField.get(sizeY-1, sizeX-1).y;

        EMGravityField.get(sizeY-1, sizeX-1).set(Ex, Ey, B/c);
    }

    public Vec2 Accel(int y, int x, float scl){
        Vec3 vF = EMGravityField.get(y, x);
        float Ex = vF.x,
              Ey = vF.y,
              B  = vF.z;

        final float h = 0.5f*dt*B, k = (1f/(1 + h*h));

        fm[0] = 1;
        fm[1] = -h;
        fm[3] = h;
        fm[4] = 1;

        Tmp.m1.set(fm).scl(k);

        Tmp.v2.set(Tmp.v1).mul(Tmp.m1);
        Tmp.v3.set(Ex, Ey).mul(Tmp.m1);
        Tmp.v4.set(Tmp.v2.y, -Tmp.v2.x).scl(0.5f*B);

        return Tmp.v2.add(Tmp.v3.add(Tmp.v4).scl(dt)).sub(Tmp.v1).scl(scl);
    }

    public void updateOther(){
        Teams.TeamData D;
        Seq<Unit> U;
        Unit u;
        float x, y, w, h;
        int ix, iy;

        for (Team t : Team.all){
            D = t.data();
            if (D == null) continue;
            U = D.units;
            for (int i = 0; i < U.size; i++){
                u = U.get(i);
                x = u.x/tilesize;
                y = u.y/tilesize;
                ix = round(x); iy = round(y);
                w = x - ix; h = y - iy;
                Tmp.v1.set(u.vel);

                if (ix < 1 || ix >= sizeX || iy < 1 || iy >= sizeY || Float.isNaN(x) || Float.isNaN(y)) continue;

                u.vel.add(Accel(iy - 1, ix - 1, (0.5f - w) * (0.5f - h)));
                u.vel.add(Accel(iy - 1, ix, (0.5f + w) * (0.5f - h)));
                u.vel.add(Accel(iy, ix - 1, (0.5f - w) * (0.5f + h)));
                u.vel.add(Accel(iy, ix, (0.5f + w) * (0.5f + h)));
            }
        }
    }

    public void initCurrent(){
        Teams.TeamData D;
        Seq<Unit> U;
        Unit u;
        float x, y, m, m2, w, h;
        int ix, iy;
        Vec2 v;
        for (Team t : Team.all) {
            D = t.data();
            if(D == null) continue;
            U = D.units;
            for(int i = 0; i < U.size; i++){
                u = U.get(i);
                x = u.x/tilesize;
                y = u.y/tilesize;
                ix = round(x); iy = round(y);
                if (ix < 1 || ix >= sizeX || iy < 1 || iy >= sizeY || Float.isNaN(x) || Float.isNaN(y)) continue;
                u.vel.limit(maxSpeed);
                v = u.vel;

                m = G * u.mass()/ tilesize / tilesize / Mathf.sqrt(1f - v.len2()/csquared);
                w = x - ix; h = y - iy;

                m2 = m * (0.5f-w)*(0.5f-h);
                CurrentField.get(iy-1, ix-1).add(csquared * m2, v.x * m2, v.y * m2);
                m2 = m * (0.5f+w)*(0.5f-h);
                CurrentField.get(iy-1, ix).add(csquared * m2, v.x * m2, v.y * m2);
                m2 = m * (0.5f-w)*(0.5f+h);
                CurrentField.get(iy, ix-1).add(csquared * m2, v.x * m2, v.y * m2);
                m2 = m * (0.5f+w)*(0.5f+h);
                CurrentField.get(iy, ix).add(csquared * m2, v.x * m2, v.y * m2);
            }
        }
    }

    public void ADI(){
        BufferField.laplacian(PotentialField, vPotentialField).sub(CurrentField).scl(0.5f*dt);
        BufferField2.laplacianY(vPotentialField).scl(0.25f*dt*dt);
        vPotentialField.add(BufferField2).add(BufferField);

        for(int i=0; i<sizeY; i++){
            tridiagonalSolver(Direction.x, i, sizeX);
        }

        BufferField2.laplacianX(vPotentialField).scl(0.25f*dt*dt);
        vPotentialField.add(BufferField2);
        for(int i=0; i<sizeX; i++){
            tridiagonalSolver(Direction.y, i, sizeY);
        }
    }

    public void tridiagonalSolver(final Direction d, int pos, final int maxSize){
        switch(d){
            case x -> tridResult[0].set(vPotentialField.get(pos, 0)).scl(b_inv);
            case y -> tridResult[0].set(vPotentialField.get(0, pos)).scl(b_inv);
        }

        for(int i=1; i<maxSize; i++){
            switch(d){
                case x -> tridResult[i].set(tridResult[i-1]).scl(-ac)
                        .add(vPotentialField.get(pos, i)).scl(tridParam2[i]);
                case y -> tridResult[i].set(tridResult[i-1]).scl(-ac)
                        .add(vPotentialField.get(i, pos)).scl(tridParam2[i]);
            }
        }

        switch(d){
            case x -> tridResult[maxSize-1].set(tridResult[maxSize-2]).scl(-2*ac)
                    .add(vPotentialField.get(pos, maxSize-1)).scl(tridParam2[maxSize-1]);
            case y -> tridResult[maxSize-1].set(tridResult[maxSize-2]).scl(-2*ac)
                    .add(vPotentialField.get(maxSize-1, pos)).scl(tridParam2[maxSize-1]);
        }

        for(int i=maxSize-2; i>=0; i--){
            tridResult[i].sub(tridResult[i+1], tridParam1[i]);
        }

        for(int i=0; i<maxSize; i++){
            switch(d){
                case x: vBuffer.get(pos, i).set(tridResult[i]);
                case y: vBuffer.get(i, pos).set(tridResult[i]);
            }
        }
    }

    public enum Direction{
        x, y;
    }

    public void updateFields(){
        vBuffer.set(vPotentialField);
        ADI();
        PotentialField.addmul(vBuffer.add(vPotentialField), 0.5f*dt);

        //Keep this one if the new one fail to work
        /*//Yoshida Integrator

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.addmul(vPotentialField, c1), dt);
        BufferField.laplacian(PotentialField, vPotentialField);
        vPotentialField.addmul(BufferField.sub(CurrentField), d1*dt);

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.addmul(vPotentialField, c2), dt);
        BufferField.laplacian(PotentialField, vPotentialField);
        vPotentialField.addmul(BufferField.sub(CurrentField), d2*dt);

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.addmul(vPotentialField, c3), dt);
        BufferField.laplacian(PotentialField, vPotentialField);
        vPotentialField.addmul(BufferField.sub(CurrentField), d3*dt);

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.addmul(vPotentialField, c4), dt);

        */

        EMGravity();

        PotentialField.normalize();
        debug = debug.concat(sizeX + ", " + sizeY);
    }

    @Override
    public int color(int i, int j){
        float s = PotentialField.get(i, j).x;
        s /= Mathf.sqrt(drawScale*drawScale + s*s);

        int color = 0xff, col;

        col = (int)(255f*Math.max(s, 0));
        color |= col << 24;

        col = (int)(255f*Math.max(-s, 0));
        color |= col << 8;
        return color;
    }

    @Override
    public void write(DataOutput write) throws IOException{
        if (PotentialField == null || vPotentialField == null) return;
        PotentialField.eachTile(v -> {
            try{
                write.writeFloat(v.x);
                write.writeFloat(v.y);
                write.writeFloat(v.z);
            }catch(IOException e){
                throw new RuntimeException(e);
            }
        });
        vPotentialField.eachTile(v -> {
            try{
                write.writeFloat(v.x);
                write.writeFloat(v.y);
                write.writeFloat(v.z);
            }catch(IOException e){
                throw new RuntimeException(e);
            }
        });
    }

    @Override
    public void read(DataInput read) throws IOException{
        sizeX = GEMgrid.sizeX = world.width();
        sizeY = GEMgrid.sizeY = world.height();

        PotentialField = new GEMgrid();
        vPotentialField = new GEMgrid();

        PotentialField.eachTile(v -> {
            try{
                v.x = read.readFloat();
                v.y = read.readFloat();
                v.z = read.readFloat();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
        vPotentialField.eachTile(v -> {
            try{
                v.x = read.readFloat();
                v.y = read.readFloat();
                v.z = read.readFloat();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
    }

    @Override
    public boolean shouldWrite(){
        return !state.rules.editor;
    }

    public static class GEMgrid{
        private final Vec3[] grid;
        private final Vec3 buffer = new Vec3();
        public static int sizeY, sizeX;

        public GEMgrid() {
            int totalSize = sizeX * sizeY;
            grid = new Vec3[totalSize];
            for(int i = 0; i < totalSize; i++){
                grid[i] = new Vec3();
            }
        }

        public final Vec3 get(int y, int x){
            return grid[x + y * sizeX];
        }

        public GEMgrid clear() {
            eachTile(v -> v.set(Vec3.Zero));
            return this;
        }

        /**
         * Note that the GEMgrid object that called this method will have its value be completely replaced.
         *
         * @param grid:  the grid that will be the values for the laplacian.
         * @param Bgrid: grid that will be used to impose the Neumann Boundary Condition: ∂(grid)/∂n = -(Bgrid)/c.
         */
        public GEMgrid laplacian(GEMgrid grid, GEMgrid Bgrid){
            for (int i = 1; i < sizeY - 1; i++) {
                for (int j = 1; j < sizeX - 1; j++){
                    get(i, j).set(grid.get(i, j)).scl(-4f)
                            .add(grid.get(i, j - 1))
                            .add(grid.get(i, j + 1))
                            .add(grid.get(i - 1, j))
                            .add(grid.get(i + 1, j))
                            .scl(csquared);
                }
            }

            for (int i = 1; i < sizeY - 1; i++) { //X-edge
                buffer.set(Bgrid.get(i, 0)).scl(-2f / c / stabilizer);
                get(i, 0).set(grid.get(i, 0)).scl(-4f)
                        .add(grid.get(i, 1))
                        .add(grid.get(i, 1))
                        .add(grid.get(i - 1, 0))
                        .add(grid.get(i + 1, 0))
                        .add(buffer).scl(csquared);

                buffer.set(Bgrid.get(i, sizeX - 1)).scl(-2f / c / stabilizer);
                get(i, sizeX - 1).set(grid.get(i, sizeX - 1)).scl(-4f)
                        .add(grid.get(i, sizeX - 2))
                        .add(grid.get(i, sizeX - 2))
                        .add(grid.get(i - 1, sizeX - 1))
                        .add(grid.get(i + 1, sizeX - 1))
                        .add(buffer).scl(csquared);
            }


            for (int j = 1; j < sizeX - 1; j++) { //Y-edge
                buffer.set(Bgrid.get(0, j)).scl(-2f / c / stabilizer);
                get(0, j).set(grid.get(0, j)).scl(-4f)
                        .add(grid.get(1, j))
                        .add(grid.get(1, j))
                        .add(grid.get(0, j - 1))
                        .add(grid.get(0, j + 1))
                        .add(buffer).scl(csquared);

                buffer.set(Bgrid.get(sizeY - 1, j)).scl(-2f / c / stabilizer);
                get(sizeY - 1, j).set(grid.get(sizeY - 1, j)).scl(-4f)
                        .add(grid.get(sizeY - 2, j))
                        .add(grid.get(sizeY - 2, j))
                        .add(grid.get(sizeY - 1, j - 1))
                        .add(grid.get(sizeY - 1, j + 1))
                        .add(buffer).scl(csquared);
            }


            //Bottom-left
            get(0, 0).set(Bgrid.get(0, 0)).scl(-1f / c / stabilizer).sub(grid.get(0, 0)).scl(2f)
                    .add(grid.get(0, 1)).add(grid.get(1, 0)).scl(2f * csquared);
            //Bottom-right
            get(0, sizeX - 1).set(Bgrid.get(0, sizeX - 1)).scl(-1f / c / stabilizer).sub(grid.get(0, sizeX - 1)).scl(2f)
                    .add(grid.get(0, sizeX - 2)).add(grid.get(1, sizeX - 1)).scl(2f * csquared);
            //Top-left
            get(sizeY - 1, 0).set(Bgrid.get(sizeY - 1, 0)).scl(-1f / c / stabilizer).sub(grid.get(sizeY - 1, 0)).scl(2f)
                    .add(grid.get(sizeY - 1, 1)).add(grid.get(sizeY - 2, 0)).scl(2f * csquared);
            //Top-right
            get(sizeY - 1, sizeX - 1).set(Bgrid.get(sizeY - 1, sizeX - 1)).scl(-1f / c / stabilizer).sub(grid.get(sizeY - 1, sizeX - 1)).scl(2f)
                    .add(grid.get(sizeY - 1, sizeX - 2)).add(grid.get(sizeY - 2, sizeX - 1)).scl(2f * csquared);

            return this;
        }

        public GEMgrid laplacianX(GEMgrid grid){
            for(int i=0; i<sizeY; i++){
                for(int j=1; j<sizeX-1; j++){
                    get(i, j).set(grid.get(i, j)).scl(-2f)
                            .add(grid.get(i, j-1))
                            .add(grid.get(i, j+1)).scl(csquared);
                }
                get(i,0).set(grid.get(i, 1)).sub(grid.get(i, 0)).scl(2f*csquared);
                get(i,sizeX-1).set(grid.get(i, sizeX-2)).sub(grid.get(i, sizeX-1)).scl(2f*csquared);
            }

            return this;
        }

        public GEMgrid laplacianY(GEMgrid grid){
            for(int j=0; j<sizeX; j++){
                for(int i=1; i<sizeY-1; i++){
                    get(i, j).set(grid.get(i, j)).scl(-2f)
                            .add(grid.get(i-1, j))
                            .add(grid.get(i+1, j)).scl(csquared);
                }
                get(0, j).set(grid.get(1, j)).sub(grid.get( 0, j)).scl(2f*csquared);
                get(sizeY-1, j).set(grid.get(sizeY-2, j)).sub(grid.get(sizeY-1, j)).scl(2f*csquared);
            }

            return this;
        }

        public GEMgrid set(GEMgrid grid){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).set(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid add(GEMgrid grid){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).add(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid sub(GEMgrid grid){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).sub(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid sub(Vec3 v){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).sub(v);
                }
            }
            return this;
        }

        public GEMgrid normalize(){
            buffer.set(Vec3.Zero);
            final float scale = 1f/grid.length;
            for(int i=1; i<sizeY-1; i++){
                buffer.add(get(i, 0))
                      .add(get(0, i))
                      .add(get(sizeY-1, i))
                      .add(get(i, sizeX-1));
            }

            buffer.add(get(0, 0))
                    .add(get(0, sizeX-1))
                    .add(get(sizeY-1, 0))
                    .add(get(sizeY-1, sizeX-1));

            buffer.scl(0.5f/((sizeX+sizeY-2)));

            sub(buffer);
            return this;
        }

        public GEMgrid scl(float s) {
            eachTile(v -> v.scl(s));
            return this;
        }

        public GEMgrid addmul(GEMgrid grid, float s){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).add(grid.get(i, j), s);
                }
            }
            return this;
        }

        public GEMgrid eachTile(Cons<Vec3> cons){
            for (Vec3 v : grid) {
                cons.get(v);
            }
            return this;
        }
    }
}
